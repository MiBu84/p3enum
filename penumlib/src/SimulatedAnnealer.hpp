/*
 * SimulatedAnnealer.hpp
 *
 *  Created on: 07.03.2019
 *      Author: mburger
 */

#ifndef SRC_SIMULATEDANNEALER_HPP_
#define SRC_SIMULATEDANNEALER_HPP_

#include <omp.h>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/progress.hpp>
#include <NTL/mat_ZZ.h>
#include <vector>
#include <iomanip>
#include <ctime>


#include "BKZBenchmarker.hpp"
#include "AnnealingSolution.hpp"
#include "Configurator.hpp"
#include "pruningfunc.hpp"

using namespace std;
using namespace NTL;

template <class FT>
class SimulatedAnnealer {
public:
	SimulatedAnnealer(const AnnealInfo<FT>& ainfo) {
		this->_ainfo = ainfo;
	}

	SimulatedAnnealer() {
		this->_ainfo = AnnealInfo<FT>();
	}

	~SimulatedAnnealer() {

	}

	int writeTempResultToFile(std::string ident, AnnealingSolution<FT>* sol_per_thread,
			const int number_of_bkz_types, const int steps, const double algstart) {
		stringstream ss;
		std::string fname = "anneal-log-";
		fname.append(ident);
		fname.append(".txt");

		ofstream myfile;
		myfile.open (fname.c_str(), ios::app);


		FT mincost = FT(std::numeric_limits<long double>::max());
		int idx = -1;
		for(int i=0; i<number_of_bkz_types;i++) {
			if(sol_per_thread[i].getCost() < mincost) {
				mincost = sol_per_thread[i].getCost();
				idx = i;
			}
		}

		ss << "Time for Annealing: " << omp_get_wtime() - algstart << " s." << endl;
		ss << "All solutions: " << endl;
		for(int i=0; i<number_of_bkz_types;i++) {
			sol_per_thread[i].printSolutionToSStream(ss);
		}

		ss << "Found final best Solution of thread " << idx <<": " << endl;
		sol_per_thread[idx].printSolutionToSStream(ss);
		ss << "Final costs: " << mincost << endl;
		ss << "After " << steps << " steps." << endl;

		myfile << ss.str();
		myfile.close();

		return 0;
	}

	void anneal(const mat_ZZ& Bi) {
		BKZBenchmarker<FT> benchi = BKZBenchmarker<FT>(_ainfo);
		mat_ZZ B = Bi;
		int dim = B.NumCols();
		vector<double> cost_prot;

		// Start with an EvenLinear function
		MBVec<double> loc_prunfunc = MBVec<double>(dim);
		for(int i=0; i<=dim-1; i++) {
			int multiplier = (i+2) / 2;
			loc_prunfunc[dim-1-i] = 1.0 * double(multiplier)/double(dim/2);
		}

		vector<FT> circ;
		for(int i=0; i<dim; i++) {
			FT tmp = ((FT)(loc_prunfunc.data()[dim-i-1]));
			circ.push_back(tmp);
			//cout << tmp << endl;
		}

		BetaConfig binfo {Configurator::getInstance().ann_prebeta_start,
			Configurator::getInstance().ann_prebeta_end,
			Configurator::getInstance().ann_beta_start,
			Configurator::getInstance().ann_beta_end};
		benchi.initialize(B, binfo);

		int numthreads = _ainfo._number_of_annealing_threads;
		// If there are more threads than BKZ-Configs, some threads will do the same Conig
		const unsigned int number_of_bkz_types = std::max(benchi.size(), (unsigned int)numthreads);

		AnnealingSolution<FT>* sol_per_thread = new AnnealingSolution<FT>[number_of_bkz_types];
		vector<double >T_init;
		T_init.reserve(number_of_bkz_types); T_init.resize(number_of_bkz_types);

		const long long numits = Configurator::getInstance().ann_iterations;
		const int init_sample_size = 100;

		double T_target = Configurator::getInstance().ann_target_temp;
		double cooling_rate = Configurator::getInstance().ann_cooling_rate;
        
        cout << "Target temp: " << T_target 
            << ", cooling rate: " << cooling_rate
            << ", iterations: " << numits << endl;
        
		long numcalcs;

		// Create a filename for the output file
		auto t = std::time(nullptr);
		auto tm = *std::localtime(&t);
		std::stringstream ss;
		ss <<  std::put_time(&tm, "%d-%m-%Y--%H-%M-%S");
		std::string starttimestr = ss.str();
		double algstart=omp_get_wtime();

		unsigned long long steps = 0;

		boost::progress_display show_progress(1);

#pragma omp parallel num_threads (numthreads)
{
		int tid = omp_get_thread_num();
		AnnealingSolution<FT> starting_sol = AnnealingSolution<FT>(&benchi, circ, _ainfo);
		AnnealingSolution<FT> act_sol;
		double t_temp_init_start=0;
		double t_temp_init_end=0;

#pragma omp for
		for(unsigned int ci = 0; ci < number_of_bkz_types; ci++) {

			// Each type is overwritten by the file
#pragma omp critical
			if(Configurator::getInstance().ann_old_result.length() > 0) {
				if(starting_sol.readSolutionFromFile(Configurator::getInstance().ann_old_result, dim, ci)!=0) {
					starting_sol = AnnealingSolution<FT>(&benchi, circ, _ainfo);
				}
			}

			// Every thread should start from the same point in a first trial
			act_sol = starting_sol;

			// Distribute the BKZ-Config around initial starting depending on thread's iteration
			act_sol.setToNumberedBetaPair(ci);
			// From here on temperature can be calculated
			t_temp_init_start = omp_get_wtime();
			FT calc_start_temp = calculateStartingTemperature(act_sol, init_sample_size);
			t_temp_init_end = omp_get_wtime();

			if(Configurator::getInstance().ann_target_temp >= 0.0) {
				T_init[ci] = 1.0;
			}
			else if (Configurator::getInstance().ann_target_temp == 2.0) {
				T_init[ci] = calc_start_temp;
			}
			else {
				T_init[ci] = Configurator::getInstance().ann_target_temp;
			}

			//cout << "T_init[" << ci << "]" << endl
			//		<< std::flush;
		}

		// Find the largest temp since it determines the highest runtime
		if(tid==0) {
			auto it = max_element(std::begin(T_init), std::end(T_init)); // c++11
			numcalcs = (log(T_target / *it)/log(cooling_rate)) * numits;

			cout << "Annealing will do " << numcalcs << " modifications per thread in estimated "
					<< (t_temp_init_end-t_temp_init_start) / (double)init_sample_size * numcalcs << " seconds.";

			show_progress.restart(numcalcs);
		}

		// Run over all considered beta-configurations
#pragma omp for
		for(unsigned int ci = 0; ci < number_of_bkz_types; ci++) {
			double T = T_init[ci];
			act_sol = starting_sol;
			act_sol.setToNumberedBetaPair(ci);

			/*cout << "Thread " << tid << " works on " << "(" << act_sol.getBetaPair().first
					<< "," << act_sol.getBetaPair().second << ")"
					<< " and ci: " << ci
					<< endl << std::flush;*/

			AnnealingSolution<FT> act_sol_back = act_sol;
			AnnealingSolution<FT> min_sol_thread = act_sol; // absolute minimal solution of thread
			sol_per_thread[ci] = min_sol_thread;

			// For decision if worse solution is accepted
			std::mt19937 rng;
			rng.seed(std::random_device()());
			std::uniform_real_distribution<double> distreal(0.0, 1.0);
			int breakcount = 0;

			long long mysteps = 0;
			long long steps_with_positive_change = 0;
			long long steps_with_negative_change = 0;

			while (T > T_target && breakcount < 1000 * numits) {
				for (int i = 0; i < numits; i++) {
					mysteps++;
					// This version of anneal always keeps its beta-config
					act_sol.modifyToNeighborSolution(false);

					// Accept solution since it is better
					if (act_sol.getCost() <= act_sol_back.getCost()) {
						act_sol_back = act_sol;
						steps_with_positive_change++;

						// No absolut minimum for thread
						if (act_sol.getCost() < min_sol_thread.getCost()) {
							min_sol_thread = act_sol;
							breakcount=0;
						}

						if(tid==0)
						{
							//cost_prot.push_back(act_sol.getCost());
						}
					}

					// Randomize if the worse solution is accepted
					else if (act_sol.getCost() > act_sol_back.getCost()) {
						FT prob = exp(-(act_sol.getCost() - act_sol_back.getCost()) / T);
						FT accept_prob = FT(distreal(rng));

						/*cout << "Prob: " << prob << " and accept prob: " << accept_prob
								<< "@ (T: " << T << ") , "
								<< "(" << act_sol.getCost() << "), "
								<< "(" << act_sol_back.getCost() << ")"
								<< endl;*/

						if (accept_prob < prob) {
							steps_with_negative_change++;
							act_sol_back = act_sol;
							breakcount=0;
							if(tid==0)
							{
								//cost_prot.push_back(act_sol.getCost());
							}
						}
						else {
							breakcount++;
							act_sol = act_sol_back;
						}
					}
					if(tid==0)
					{
						++show_progress;
						steps++;
					}
				}
				T = T * cooling_rate;
			}
			sol_per_thread[tid] = min_sol_thread;
			cout << "This run with: " << mysteps << "(Positive: "
					<< steps_with_positive_change << " / Negative: "
					<< steps_with_negative_change << ")" << endl;
		}

} // End parallel

		FT mincost = FT(std::numeric_limits<long double>::max());
		int idx = -1;
		for(unsigned int i=0; i<number_of_bkz_types;i++) {
			if(sol_per_thread[i].getCost() < mincost) {
				mincost = sol_per_thread[i].getCost();
				idx = i;
			}
		}

		cout << "Time for Annealing: " << omp_get_wtime() - algstart << " s." << endl;
		cout << "All solutions: " << endl;
		for(unsigned int i=0; i<number_of_bkz_types;i++) {
			sol_per_thread[i].printSolution();
		}


		cout << "Found final best Solution of thread " << idx <<": " << endl;
		sol_per_thread[idx].printSolution();
		cout << "Final costs: " << mincost << endl;
		cout << "After " << steps << " steps." << endl;

		writeTempResultToFile(starttimestr, sol_per_thread, number_of_bkz_types, steps, algstart);


		/*ofstream myfile;
		myfile.open ("Lanbaer.txt", ios::app);

		int i=0;
		for(auto it = cost_prot.begin(); it != cost_prot.end();++it) {
			myfile << i << "\t" << *it << "\n";
			i++;
		}

		myfile.close();*/

		delete[] sol_per_thread;
		exit(222);
	}






						//algend=omp_get_wtime();
//					}

					/*if(algend - algstart > 300) {
						sol_per_thread[tid] = min_sol_thread;
#pragma omp barrier
#pragma omp master
						writeTempResultToFile(starttimestr, sol_per_thread, numthreads, steps);
#pragma omp barrier
#pragma omp master
						algstart=algend;
					}*/
//				}

//			}

//		}



//	}

private:
	FT calculateStartingTemperature(AnnealingSolution<FT>& input_sol, int sample_size=100) {
		vector<AnnealingSolution<FT>> sol_array;
		sol_array.reserve(sample_size);

		// Generate random neighboring solutions
		for(int i=0; i<sample_size; i++) {
			AnnealingSolution<FT> vary = input_sol;
			vary.modifyToNeighborSolution(false);
			sol_array.push_back(vary);
		}

		// Find biggest distance between all samples
		FT max_dist = FT(0.0);
		FT act_dist = FT(0.0);
		for(int i=0; i<sample_size; i++) {
			for(int j=0; j<sample_size; j++) {
				act_dist = abs(sol_array[i].getCost() - sol_array[j].getCost());

				if(act_dist > max_dist)
					max_dist=act_dist;
			}
		}

		FT start_temp = -max_dist / (log(0.8));
		cout << "Calculated starting temperature of " << start_temp << " (max_dist: " << max_dist << ")" << endl;
		return start_temp;
	}

	AnnealInfo<FT> _ainfo;
};

#endif /* SRC_SIMULATEDANNEALER_HPP_ */
