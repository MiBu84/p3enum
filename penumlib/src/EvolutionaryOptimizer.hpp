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

template <class FT>
class EvolutionaryOptimizer {
public:
	EvolutionaryOptimizer(const AnnealInfo<FT>& ainfo) {
		this->_ainfo = ainfo;
		_threads=ainfo._number_of_annealing_threads;

		for(int i=0; i<_threads; i++)
			_populations.push_back(EvolutionPopulation<FT>(200));

	}

	EvolutionaryOptimizer() {
		this->_ainfo = AnnealInfo<FT>();
		_populations.push_back(EvolutionPopulation<FT>(200));
		_threads=1;
	}

	EvolutionaryOptimizer(EvolutionaryOptimizer& other) {
		this->_ainfo = other._ainfo;
		this->_threads = other._threads;
		for(int i=0; i<_threads; i++)
			this->_populations[i] = other._populations[i];
	}

	EvolutionaryOptimizer(const EvolutionaryOptimizer& other) {
		this->_ainfo = other._ainfo;
		this->_threads = other._threads;
		for(int i=0; i<_threads; i++)
			this->_populations[i] = other._populations[i];
	}

	EvolutionaryOptimizer& operator= (const EvolutionaryOptimizer & other) {
		this->_ainfo = other._ainfo;
		this->_threads = other._threads;
		for(int i=0; i<_threads; i++)
			this->_populations[i] = other._populations[i];

		return *this;
	}


	void doEvolution(const mat_ZZ& Bi, MBVec<double>& pruningfuction) {
		cout << "Starting Evolution" << endl;

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

		const int generations = 400;
		boost::progress_display show_progress(generations);

#pragma omp parallel num_threads (numthreads)
{
		int tid = omp_get_thread_num();
		//cout << "My tid: " << tid << endl;

#pragma omp for
		for(unsigned int ci = 0; ci < number_of_bkz_types; ci++) {
			EvolutionarySolution<FT> starting_sol = EvolutionarySolution<FT>(&benchi, circ, _ainfo);
			starting_sol.setToNumberedBetaPair(ci);

			int attempcnt = 0;
			int succcnt = 0;

			EvolutionarySolution<FT> tomodify_sol = starting_sol;

			// Generate target_size random solutions
			_populations[tid].clear();
			_populations[tid].insertIndividual(tomodify_sol);

			for(unsigned int i=1; i < _populations[tid].getTargetSize(); i++) {
				for(int ii=0; ii<10000; ii++) {
					tomodify_sol.modifyToNeighborSolution(false, false);
				}
				tomodify_sol.calculateCosts();
				_populations[tid].insertIndividual(tomodify_sol);
				tomodify_sol = starting_sol;

			}

			for (int ac =0; ac<generations; ac++) {
				FT best = _populations[tid].getMaxFitness();
				_populations[tid].nextGeneration();
				FT new_best = _populations[tid].getMaxFitness();

				if (new_best > best) {
					succcnt++;
				}
				attempcnt++;

				if(tid==0)
					++show_progress;
			}

			cout << succcnt << " out of " << attempcnt << endl;
		} // End BKZ types

#pragma omp single
		{
			for(int i=0; i<_threads; i++) {
				_populations[i].getBestIndiviual().printSolution();
			}
		}

	} // end parallel
	}




protected:
	AnnealInfo<FT> _ainfo;
	std::vector<EvolutionPopulation<FT>> _populations;
	int _threads;
};


#endif /* SRC_EVOLUTIONARYOPTIMIZER_HPP_ */
