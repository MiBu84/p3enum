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

template <class FT>
class EvolutionaryOptimizer {
public:
	EvolutionaryOptimizer(const AnnealInfo<FT>& ainfo) {
		this->_ainfo = ainfo;
		_population = EvolutionPopulation<FT>(10);
	}

	EvolutionaryOptimizer() {
		this->_ainfo = AnnealInfo<FT>();
		_population = EvolutionPopulation<FT>(10);
	}

	EvolutionaryOptimizer(EvolutionaryOptimizer& other) {
		this->_ainfo = other._ainfo;
		this->_population = other._population;
	}

	EvolutionaryOptimizer(const EvolutionaryOptimizer& other) {
		this->_ainfo = other._ainfo;
		this->_population = other._population;
	}

	EvolutionaryOptimizer& operator= (const EvolutionaryOptimizer & other) {
		this->_ainfo = other._ainfo;
		this->_population = other._population;
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

		EvolutionarySolution<FT> starting_sol = EvolutionarySolution<FT>(&benchi, circ, _ainfo);
		starting_sol.setToNumberedBetaPair(0);
		starting_sol.printCosts();

		starting_sol.modifyToNeighborSolution(false, false);
		starting_sol.calculateCosts();
		starting_sol.printCosts();

		EvolutionarySolution<FT> tomodify_sol = starting_sol;

		for(int i=0; i<10000; i++) {
			tomodify_sol.modifyToNeighborSolution(false, false);

			if(i%1000==0) {
				tomodify_sol.calculateCosts();
				_population.insertIndividual(tomodify_sol);
			}
		}
		//_population.printCostsOfIndividuals();

		//_population.printSelectionProbablities();
		//cout << "End" << endl;


		tomodify_sol.calculateCosts();
		tomodify_sol.printCosts();

		//EvolutionarySolution<FT> res_of_recombine = recombine(starting_sol, tomodify_sol);
		EvolutionarySolution<FT> res_of_recombine = crossingOver(starting_sol, tomodify_sol);
		res_of_recombine.printCosts();

		_population.printCostsOfIndividuals();
		cout << "New step " << endl;
		_population.nextGeneration();
		_population.printCostsOfIndividuals();

		/*map <int, int> mappi;
		cout << "Thrown dice: " << endl;
		for (int i=0; i<50000; i++) {
			int id = _population.selectionWithRoullette().getRandomID();

			auto it = mappi.find(id);
			if(it==mappi.end()) {
				mappi[id] = 0;
			}

			else {
				(it->second)++;
			}
		}

		for (auto it=mappi.begin(); it != mappi.end(); it++) {
			cout << it->first << " " << it->second << endl;
		}
		cout << endl;

		_population.printCostsOfIndividuals();*/
	}



private:
	/**
	 * Recombines the pruning functions of two different solutions by cutting at a legal position and
	 * appending the second part of sol2 to sol1
	 */
	EvolutionarySolution<FT> recombine(const EvolutionarySolution<FT>& sol1, const EvolutionarySolution<FT>& sol2) {
		EvolutionarySolution<FT> sol_res = EvolutionarySolution<FT>(sol1);

		// Find out where to cut
		std::random_device seeder_mod;
		std::mt19937 engine(seeder_mod());
		std::uniform_int_distribution<int> dist_dim(1, (sol1._dim-2)/2);
		int cut_dim = dist_dim(engine);

		cout << "Cut at " << cut_dim*2 << endl;

		// Now copy the corresponding positions of pruning function
		for(int i=cut_dim*2; i < sol1._dim-2; i++) {
			sol_res._prun_func[i] = sol2._prun_func[i];
		}

		// Invalid solutions have to be repaired
		for(int i=cut_dim*2; i < sol1._dim-2; i++) {
			if(sol_res._prun_func[i] > sol_res._prun_func[i+1]) {
				sol_res._prun_func[i+1] = sol_res._prun_func[i];
			}
			else {
				break;
			}
		}

		// Recalc probs and costs!
		sol_res.resetProbs();
		sol_res.calculateCosts();

		return sol_res;
	}

	EvolutionarySolution<FT> crossingOver(const EvolutionarySolution<FT>& sol1, const EvolutionarySolution<FT>& sol2) {
		EvolutionarySolution<FT> sol_res = EvolutionarySolution<FT>(sol1);

		// Find out range to crossover
		std::random_device seeder_mod;
		std::mt19937 engine(seeder_mod());

		std::uniform_int_distribution<int> dist_dim_start(0, (sol1._dim-2)/2);
		int start_dim = dist_dim_start(engine);

		std::uniform_int_distribution<int> dist_dim_end(start_dim+1, (sol1._dim-1)/2);
		int end_dim  = dist_dim_end(engine);

		if(start_dim >= end_dim || end_dim*2+1 >= sol1._dim) {
			int stop=0;
			cerr << "EvolutionarySolutionCrossing over with start " << start_dim
					<< " and end " << end_dim << endl;
			cin >> stop;
		}

		cout << "EvolutionarySolutionCrossing over with start " << start_dim
				<< " and end " << end_dim << endl;

		// Copy the values
		for(int i = start_dim*2; i <= end_dim*2+1; i++) {
			sol_res._prun_func[i] = sol2._prun_func[i];
		}

		// Invalid solutions have to be repaired
		FT low_crossval = sol2._prun_func[start_dim*2];
		FT high_crossval = sol2._prun_func[end_dim*2];

		for(int i=start_dim*2; i>=1; i--) {
			if(low_crossval < sol_res._prun_func[i-1]) {
				sol_res._prun_func[i-1] = low_crossval;
			}
			else {
				break;
			}
		}

		for(int i=end_dim*2; i <= sol_res._dim-3; i++) {
			if(high_crossval > sol_res._prun_func[i+1]) {
				sol_res._prun_func[i+1] = high_crossval;
			}
			else {
				break;
			}
		}

		// Recalc probs and costs!
		sol_res.resetProbs();
		sol_res.calculateCosts();
		return sol_res;
	}
protected:
	AnnealInfo<FT> _ainfo;
	EvolutionPopulation<FT> _population;
};


#endif /* SRC_EVOLUTIONARYOPTIMIZER_HPP_ */
