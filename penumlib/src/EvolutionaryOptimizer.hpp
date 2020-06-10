/*
 * EvolutionaryOptimizer.hpp
 *
 *  Created on: 15.10.2019
 *      Author: mburger
 */

#ifndef SRC_EVOLUTIONARYOPTIMIZER_HPP_
#define SRC_EVOLUTIONARYOPTIMIZER_HPP_

#include "EvolutionarySolution.hpp"
#include "EvolutionPopulation.hpp"
#include "BKZBenchmarker.hpp"
#include "MBVec.hpp"
#include <boost/progress.hpp>
//#include <boost/timer/progress_display.hpp>

template <class FT>
class EvolutionaryOptimizer {
public:
	EvolutionaryOptimizer(const AnnealInfo<FT>& ainfo) {
		this->_ainfo = ainfo;
		_threads=ainfo._number_of_annealing_threads;

		/*for(int i=0; i<_threads; i++)
			_populations.push_back(EvolutionPopulation<FT>(250));*/

	}

	EvolutionaryOptimizer() {
		this->_ainfo = AnnealInfo<FT>();
		/*_populations.push_back(EvolutionPopulation<FT>(250));*/
		_threads=1;
	}

	EvolutionaryOptimizer(EvolutionaryOptimizer& other) {
		this->_ainfo = other._ainfo;
		this->_threads = other._threads;
		/*for(int i=0; i<_threads; i++)
			this->_populations[i] = other._populations[i];*/
	}

	EvolutionaryOptimizer(const EvolutionaryOptimizer& other) {
		this->_ainfo = other._ainfo;
		this->_threads = other._threads;
		/*for(int i=0; i<_threads; i++)
			this->_populations[i] = other._populations[i];*/
	}

	EvolutionaryOptimizer& operator= (const EvolutionaryOptimizer & other) {
		this->_ainfo = other._ainfo;
		this->_threads = other._threads;
		/*for(int i=0; i<_threads; i++)
			this->_populations[i] = other._populations[i];*/

		return *this;
	}


	void doEvolution(const mat_ZZ& Bi, MBVec<double>& pruningfuction) {
		cout << "Starting Evolution" << endl;
		auto t = std::time(nullptr);
		auto tm = *std::localtime(&t);
		std::stringstream ss;
		ss <<  std::put_time(&tm, "%d-%m-%Y--%H-%M-%S");
		std::string starttimestr = ss.str();

		BKZBenchmarker<FT> benchi = BKZBenchmarker<FT>(_ainfo);
		mat_ZZ B = Bi;
		int dim = B.NumCols();
		vector<double> cost_prot;

		// Start with an EvenLinear function
		//Todo: Unfify both loops to one
		MBVec<double> loc_prunfunc =  pruningfuction;
		for(int i=0; i<=dim-1; i++) {
			int multiplier = (i+2) / 2;
			loc_prunfunc[dim-1-i] = 1.0 * double(multiplier)/double(dim/2);
		}

		vector<FT> circ;
		for(int i=0; i<dim; i++) {
			FT tmp = ((FT)(loc_prunfunc.data()[dim-i-1]));
			circ.push_back(tmp);
		}

		BetaConfig binfo {Configurator::getInstance().ann_prebeta_start,
			Configurator::getInstance().ann_prebeta_end,
			Configurator::getInstance().ann_beta_start,
			Configurator::getInstance().ann_beta_end};
		benchi.initialize(B, binfo);

		// From here threads have to do different work
		int numthreads = _ainfo._number_of_annealing_threads;
		const unsigned int number_of_bkz_types = std::max(benchi.size(), (unsigned int)numthreads);
		srand(time(0));

		const int generations = Configurator::getInstance().evo_generations;
		const int generations_to_process = generations*number_of_bkz_types;
		int generations_done = 0;
		const int one_percent_of_generations = (int)ceil(((double)generations_to_process*0.01));
		int next_generation_target = one_percent_of_generations;
		int print_cnt = 0;
		int stages_without_improvement = 0;
		//boost::progress_display show_progress(generations);


		double tstart = omp_get_wtime();
		EvolutionarySolution<FT> starting_sol_glob = EvolutionarySolution<FT>(&benchi, circ, _ainfo);
		starting_sol_glob.calculateCosts();

		EvolutionarySolution<FT> glob_best_sol = starting_sol_glob;
		cout << "Initial costs: " << starting_sol_glob.getCost() << endl;

#pragma omp parallel num_threads (numthreads) shared(generations_done)
{
		int tid = omp_get_thread_num();
		EvolutionPopulation<FT> mypop = EvolutionPopulation<FT>(300);
		//.push_back(EvolutionPopulation<FT>(250));
		EvolutionarySolution<FT> starting_sol_loc = EvolutionarySolution<FT>(&benchi, circ, _ainfo);

#pragma omp for
		for(unsigned int ci = 0; ci < number_of_bkz_types; ci++) {

			starting_sol_loc.setToNumberedBetaPair(ci);

			int attempcnt = 0;
			int succcnt = 0;

			EvolutionarySolution<FT> tomodify_sol = starting_sol_loc;
            
            // Warmup

			// Generate target_size random solutions
			mypop.clear();
			mypop.insertIndividual(tomodify_sol);

			for(unsigned int i=1; i < mypop.getTargetSize(); i++) {
				for(int ii=0; ii<1000; ii++) {
					tomodify_sol.modifyToNeighborSolution(false, false);
				}
				tomodify_sol.calculateCosts();
				mypop.insertIndividual(tomodify_sol);
 				tomodify_sol = starting_sol_loc;
			}
            
            #pragma omp critical 
            {
                if(glob_best_sol.getCost() > mypop.getBestIndiviual().getCost()) {
                    glob_best_sol = mypop.getBestIndiviual();
                }
            }

			for (int ac =0; ac<generations &&
			stages_without_improvement < one_percent_of_generations*6;
			ac++) {
				FT best = mypop.getMaxFitness();
				mypop.nextGeneration();
				FT new_best = mypop.getMaxFitness();

				if (new_best > best) {
					succcnt++;
				}
				attempcnt++;

#pragma omp critical
				{
					if(glob_best_sol.getCost() > mypop.getBestIndiviual().getCost()) {
						// Calculate improvement
						FT improvement = (glob_best_sol.getCost() - mypop.getBestIndiviual().getCost()) - FT(1.0);

						if(abs(improvement) < 5e-3) {
							stages_without_improvement++;
							cout << "Impr:" << improvement << endl;
						}
						else
							stages_without_improvement=0;

						glob_best_sol = mypop.getBestIndiviual();


					}
					else {
						stages_without_improvement++;
					}
				}

				// Verbosity
#pragma omp atomic
				generations_done++;

				if(tid==0) {
					if(generations_done > next_generation_target) {
						next_generation_target += one_percent_of_generations;
						print_cnt++;

						double tend_consumed = omp_get_wtime() - tstart;
						double time_per_gen = tend_consumed / generations_done;
						int generations_remaining = generations_to_process - generations_done;
						double expected_remaining_time = time_per_gen * (double)generations_remaining;

						cout << "(" << print_cnt << " \%): " << generations_done << "/" << generations_to_process
								<< ", lowest costs: " << glob_best_sol.getCost()
								<< ", remaining: " << expected_remaining_time << "s." << endl;
					}
					//++show_progress;
				}
			}

			cout << succcnt << " out of " << attempcnt << endl;
            //mypop.printStatistics();
		} // End BKZ types

#pragma omp critical
		{
			/*writeTempResultToFile(starttimestr, mypop, number_of_bkz_types, generations_done);

			FT mincost = numeric_limits<FT>::max();
			unsigned int mincost_idx = mypop.size();

			for(int i=0; i<_threads; i++) {
				if(_populations[i].getBestIndiviual().getCost() < mincost) {
					mincost = _populations[i].getBestIndiviual().getCost();
					mincost_idx = i;
				}
			}

			if (mincost_idx < _populations.size()) {

				_populations[mincost_idx].getBestIndiviual().printSolution();
				_populations[mincost_idx].printStatistics();
			}

			else {
				cerr << "Index " << mincost_idx << " for minimal solution invalid." << endl;
			}*/

		}

	} // end parallel
	cout << "Found minimal solution: " << endl;
	glob_best_sol.printSolution();

	//writeTempResultToFile(starttimestr, mypop, number_of_bkz_types, generations_done);

	cout << "Evolution's duration: " << omp_get_wtime() - tstart << " sec." << endl;
	}




protected:
	AnnealInfo<FT> _ainfo;
	//std::vector<EvolutionPopulation<FT>> _populations;
	int _threads;

	int writeTempResultToFile(std::string ident, EvolutionPopulation<FT> pops,
				const int numthreads, const int generation) {
			std::string fname = "evolution-log-";
			fname.append(ident);
			fname.append(".txt");

			ofstream myfile;
			myfile.open (fname.c_str(), ios::app);


			FT mincost = FT(std::numeric_limits<long double>::max());
			int idx = -1;
			/*for(int i=0; i<numthreads;i++) {
				if(pops.getBestIndiviual().getCost() < mincost) {
					mincost = pops[i].getBestIndiviual().getCost();
					idx = i;
				}
			}*/

			stringstream ss;

			ss << "Found final Solution of thread " << idx <<": " << endl;
			pops.getBestIndiviual().printSolutionToSStream(ss);
			ss << "Final costs: " << mincost << endl;
			ss << "In " << generation << " generations." << endl;

			myfile << ss.str();
			myfile.close();

			return 0;

		}
};


#endif /* SRC_EVOLUTIONARYOPTIMIZER_HPP_ */
