/*
 * EvolutionPopulation.hpp
 *
 *  Created on: 16.10.2019
 *      Author: mburger
 */

#ifndef SRC_EVOLUTIONPOPULATION_HPP_
#define SRC_EVOLUTIONPOPULATION_HPP_

#include <set>

#include "EvolutionarySolution.hpp"

#ifndef m_assert
#define m_assert(expr, msg) assert(( (void)(msg), (expr) ))
#endif

using namespace std;
//template <class FT>
//using individual = std::pair <int, EvolutionarySolution<FT>> ;

template <class FT>
struct classcomp {
	bool operator() (const EvolutionarySolution<FT>& lhs, const EvolutionarySolution<FT>& rhs) const
	{
		if (lhs.getFitness() > rhs.getFitness()) {
			return true;
		}

		else if (lhs.getFitness() == rhs.getFitness() && lhs.getRandomID() > rhs.getRandomID()) {
			return true;
		}
		else {
			return false;
		}
		//return lhs.getFitness() > rhs.getFitness();
		//return lhs.getRandomID() > rhs.getRandomID();
	}
};

template <class FT>
struct classcompID {
	bool operator() (const EvolutionarySolution<FT>& lhs, const EvolutionarySolution<FT>& rhs) const
	{
		return lhs.getRandomID() > rhs.getRandomID();
	}
};


template <class FT>
class EvolutionPopulation {
public:
	EvolutionPopulation(unsigned int target_size) {
		this->_population.clear();
		this->_target_size = target_size;
		//cout << "Setting target size to " << this->_target_size << endl;
		this->_best_individual = EvolutionarySolution<FT>();
		this->_max_costs = 0;
		this->no_of_mutations = this->no_of_shifts = this->no_of_recombines = this->no_of_crossing_over = this->no_of_keeping = 0;
	}

	EvolutionPopulation() {
		this->_population.clear();
		this->_target_size = 101;
		//cout << "Setting default target size to " << this->_target_size << endl;
		this->_best_individual = EvolutionarySolution<FT>();
		this->_max_costs = 0;
		this->no_of_mutations = this->no_of_shifts = this->no_of_recombines = this->no_of_crossing_over = this->no_of_keeping = 0;
	}

	~EvolutionPopulation () {
		this->no_of_mutations = this->no_of_shifts = this->no_of_recombines = this->no_of_crossing_over = this->no_of_keeping = 0;
	}

	EvolutionPopulation(EvolutionPopulation& other) {
		this->_target_size = other._target_size;
		this->_max_costs = other._max_costs;
		this->_best_individual = other._best_individual;
		this->_population = other._population;
		this->no_of_shifts = other.no_of_shifts;
		this->no_of_recombines = other.no_of_recombines;
		this->no_of_crossing_over = other.no_of_crossing_over;
		this->no_of_keeping = other.no_of_keeping;
		this->no_of_mutations = other.no_of_mutations;
	}

	EvolutionPopulation(const EvolutionPopulation& other) {
		this->_target_size = other._target_size;
		this->_max_costs = other._max_costs;
		this->_best_individual = other._best_individual;
		this->_population = other._population;
		this->no_of_shifts = other.no_of_shifts;
		this->no_of_recombines = other.no_of_recombines;
		this->no_of_crossing_over = other.no_of_crossing_over;
		this->no_of_keeping = other.no_of_keeping;
		this->no_of_mutations = other.no_of_mutations;
	}

	EvolutionPopulation& operator= (const EvolutionPopulation & other) {
		_target_size = other._target_size;
		_max_costs = other._max_costs;
		_best_individual = other._best_individual;
		_population = other._population;
		no_of_shifts = other.no_of_shifts;
		no_of_recombines = other.no_of_recombines;
		no_of_crossing_over = other.no_of_crossing_over;
		no_of_keeping = other.no_of_keeping;
		no_of_mutations = other.no_of_mutations;
		return *this;
	}

	EvolutionPopulation& operator[] (int index) {
		auto it = this->_population.begin();
		std::advance(it, index);
		return *it;
	}

	void clear() {
		this->_population.clear();
	}

	int insertIndividual (EvolutionarySolution<FT>& indi) {
		EvolutionarySolution<FT> indi_copy = indi;
		this->_population.insert(indi_copy);

		if(indi.getFitness() > this->_best_individual.getFitness()) {
			_best_individual = indi;
		}

		return this->_population.size();
	}

	EvolutionarySolution<FT> returnRandomIndividual() {

	}

	int reducePopulation () {
		return 0;
	}

	FT nextGeneration() {
		std::set<EvolutionarySolution<FT>, classcomp<FT>> next_gen;
		std::set<EvolutionarySolution<FT>, classcomp<FT>> parents;

		unsigned int parent_size = _target_size;

		// Choose those that should breed
		int rank=0;
		while (parents.size() != parent_size) {
			// Make copy with new ID
			//EvolutionarySolution<FT> tmp = selectionWithRoullette();
			int access = (int)floor(rank/2);
			//cout <<  access << endl;
			EvolutionarySolution<FT> tmp = selectionByRank(access);
			parents.insert( tmp );
			rank++;
		}

		/*for(auto it=parents.begin(); it != parents.end();++it) {
			cout << it->getFitness() << endl;
		}
		cout << "End parents" << endl;*/

		while (next_gen.size() != _target_size) {
		//for(unsigned int elemcnt=0; elemcnt < _target_size; elemcnt++) {

			// The next generation should be of same size
			int choser = rand() % 5;
			EvolutionarySolution<FT> new_gen_child;

			// Case 0: A random parent will also be in the next generation
			if(choser==0) {
				auto it = parents.begin();
				unsigned int rnd = rand() % (parent_size - 1);

				m_assert(rnd < parent_size  - 1, "rnd must be equal larger than parent_size");

				std::advance(it, rnd);
				new_gen_child = *it;
				this->no_of_keeping++;
			}

			// Case 1: A parent is mutated in a random amount
			else if(choser==1) {
				auto it = parents.begin();
				unsigned int rnd = (rand() % (parent_size  - 1));

				m_assert(rnd < parent_size  - 1, "rnd must be equal larger than parent_size");

				std::advance(it, rnd);
				new_gen_child = *it;

				int maxval = rand()%10+1;
				for(int i = 0 ; i < maxval; i++) {
					new_gen_child.modifyToNeighborSolution(false, false);
				}
				new_gen_child.calculateCosts();
				this->no_of_mutations++;
			}

			// Case 2: Two random solutions do CrossingOver
			else if(choser==2) {
				auto it1 = parents.begin();
				auto it2 = parents.begin();

				unsigned int rnd1 = rand() % (parent_size  - 1);
				unsigned int rnd2 = rand() % (parent_size  - 1);

				m_assert(rnd1 < parent_size  - 1, "rnd1 must be equal larger than parent_size");
				m_assert(rnd2 < parent_size  - 1, "rnd2 must be equal larger than parent_size");

				std::advance(it1, rnd1);
				std::advance(it2, rnd2);

				new_gen_child =  crossingOver(*it1, *(it2));
				this->no_of_crossing_over++;
			}

			// Case 3: Two random solutions do Recombination
			else if(choser==3) {
				auto it1 = parents.begin();
				auto it2 = parents.begin();

				unsigned int rnd1 = rand() % (parent_size  - 1);
				unsigned int rnd2 = rand() % (parent_size  - 1);

				m_assert(rnd1 < parent_size  - 1, "rnd1 must be equal larger than parent_size");
				m_assert(rnd2 < parent_size  - 1, "rnd2 must be equal larger than parent_size");

				std::advance(it1, rnd1);
				std::advance(it2, rnd2);

				new_gen_child =  recombine(*it1, *(it2));
				this->no_of_recombines++;
			}

			// Case 4: Shift one solution
			else if(choser==4) {
				auto it = parents.begin();
				unsigned int rnd = (rand() % (parent_size  - 1));

				m_assert(rnd < parent_size  - 1, "rnd must be equal larger than parent_size");

				std::advance(it, rnd);
				new_gen_child = shift(*it);
				this->no_of_shifts++;
			}

			if(new_gen_child.getFitness() > this->_best_individual.getFitness()) {
				this->_best_individual = new_gen_child;
			}


			next_gen.insert(new_gen_child);
			_population = next_gen;
		}

		if(_population.size() != _target_size) {
			int stop = 0;
			cerr << "Pop size is " << _population.size() << " but should be "
					<< _target_size << endl;
			cin >> stop;
		}

		return 0;
	}

	void printFitnessOfIndividuals() {
		auto iti = _population.begin();

		for(;iti != _population.end(); ++iti) {
			cout << "ID(" << iti->getRandomID() << ")"
			<< iti->getFitness() << endl;
		}

		cout << "Best individual: " << this->_best_individual.getFitness() << endl;
		cout << "End of costs." << endl;


		/*iti = _population.begin();
		cout << "Highest val= " << iti->getFitness() << endl;
		iti = _population.end(); iti--;
		cout << "Lowest val= " << iti->getFitness() << endl;*/
	}

	void printSelectionProbablities() {
		// Get the worst element
		FT sum_fitness = FT(0);

		auto iti = _population.begin();
		for(; iti != _population.end(); ++iti) {
			sum_fitness+= iti->getFitness();
		}

		FT sum_probs = FT(0);
		iti = _population.begin();
		for(; iti != _population.end(); ++iti) {
			sum_probs+= iti->getFitness() / sum_fitness;
			cout << "ID(" << iti->getRandomID() << "): ";
			cout << iti->getFitness() / sum_fitness << endl;
		}
		cout << "End of probs." << endl;
	}

	void printStatistics() {
		cout << "Population stastics: " << endl
			<< "No of Mutations: " << no_of_mutations << endl
			<< "No of Shifts: " << no_of_shifts << endl
			<< "No of Recombines: " << no_of_recombines << endl
			<< "No_of_CrossingOvers: " << no_of_crossing_over << endl
			<< "No_of_Survives: " << no_of_keeping << endl;
	}

	EvolutionarySolution<FT> selectionWithRoullette(const std::set<EvolutionarySolution<FT>> pop=NULL) {
		// Set a random id
		std::random_device seeder;
		std::mt19937 engine(seeder());
		std::uniform_real_distribution<FT> dist(0.0, 1.0);

		FT sum_fitness = FT(0);
		auto iti = _population.begin();
		for(; iti != _population.end(); ++iti) {
			sum_fitness+= iti->getFitness();
		}

		FT prob = dist(engine);

		FT sum_probs = FT(0);

		if(pop==NULL) {
			auto it = _population.begin();
			for(; it != _population.end(); ++it) {
				if(prob <= sum_probs) {
					return *it;
				}
				else {
					sum_probs+= (it->getFitness() / sum_fitness);
					if(prob <= sum_probs) {
						return *it;
					}
				}
			}
			cerr << "Err: " << prob << " / " << sum_probs  << endl;
			return *it;
		}

		else {
			auto it = pop.begin();
			for(; it != pop.end(); ++it) {
				if(prob <= sum_probs) {
					return *it;
				}
				else {
					sum_probs+= (it->getFitness() / sum_fitness);
					if(prob <= sum_probs) {
						return *it;
					}
				}
			}
			cerr << "Err: " << prob << " / " << sum_probs  << endl;
			return *it;
		}
	}

	EvolutionarySolution<FT> selectionByRank(int rank) {
		auto it = this->_population.begin();
		std::advance(it, rank);
		return *it;
	}

	unsigned int getTargetSize() {
		return this->_target_size;
	}

	FT getMaxFitness() {
		return this->_best_individual.getFitness();
	}

	EvolutionarySolution<FT>& getBestIndiviual() {
		return this->_best_individual;
	}

protected:
	std::set<EvolutionarySolution<FT>, classcomp<FT>> _population; /**/
	EvolutionarySolution<FT> _best_individual;
	double _max_costs;
	unsigned int _target_size;

private:
	int no_of_shifts, no_of_recombines, no_of_crossing_over, no_of_keeping, no_of_mutations;

	EvolutionarySolution<FT> shift(const EvolutionarySolution<FT>& sol1) {
		EvolutionarySolution<FT> sol_res = EvolutionarySolution<FT>(sol1);
		int shiftint = (rand()%201) - 100; // between -100 and 100
		FT shiftfloat = (FT)shiftint / FT(1000.0); // shift between -0.01 und 0.01

		for(int i=0; i < sol1._dim; i++) {
			// Prevent shift to negative values
			sol_res._prun_func[i] = std::min<FT>(std::max<FT> (sol_res._prun_func[i] + shiftfloat, FT(1e-3)), FT(1.0));
		}

		sol_res.resetProbs();
		sol_res.calculateCosts();
		return sol_res;
	}

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

		// Now copy the corresponding positions of pruning function
		for(int i=cut_dim*2; i < sol1._dim-2; i++) {
			sol_res._prun_func[i] = sol2._prun_func[i];
		}

		// Invalid solutions have to be repaired
		//ToDO: Copy the highest value from before towards the higher entries, until function is valid
		for(int i=cut_dim*2; i > 0; i--) {
			if(sol_res._prun_func[i] < sol_res._prun_func[i-1]) {
				sol_res._prun_func[i-1] = sol_res._prun_func[i];
				// Slight increment because of numerical stability
				/*if(i%2!=0)
					sol_res._prun_func[i+1] += sol_res._prun_func[i+1] * 1e-8;*/
			}
			else {
				break;
			}
		}

		// Check feasiability
		/*for(long unsigned int i=0; i<sol_res._prun_func.size()-1; i++) {
				if(sol_res._prun_func[i] > sol_res._prun_func[i+1]) {
					cout << "Prunfunc invalid! :" << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
					<< i << " (" << sol_res._prun_func[i] << " / " << sol_res._prun_func[i+1] << ")" << endl << std::flush;


					cout << "EvolutionarySolutionRecombine with cut at " << cut_dim*2
							<< " with ID(" << sol1.getRandomID() << ") and ID("
							<< sol2.getRandomID() << ")"
							<< endl << std::flush;

					// Copy the values
					cout << "Start writing: " << endl;
					for(int i = cut_dim*2; i < sol1._dim-2; i++) {
						sol_res._prun_func[i] = sol2._prun_func[i];
						cout << "[" << i << "]: " << sol2._prun_func[i] << endl;
					}
					cout << "End writing" << endl;

					sol_res.printSolution(true);

					cout << endl << "Above: ";
					// Check above
					FT high_crossval = sol_res._prun_func[cut_dim*2+1];
					for(int i=cut_dim*2+1; i < sol1._dim-2; i++) {
						cout << high_crossval << "/" << sol_res._prun_func[i+1] << " ";
						if(high_crossval > sol_res._prun_func[i+1]) {
							cout << i << endl;
						}
						else {
							break;
						}
					}
					cout << endl;


					exit(-1000);
				}
		}*/

		// Recalc probs and costs!
		sol_res.resetProbs();
		sol_res.calculateCosts();

		return sol_res;
	}

	EvolutionarySolution<FT> crossingOver(const EvolutionarySolution<FT>& sol1, const EvolutionarySolution<FT>& sol2) {
		EvolutionarySolution<FT> sol_res = EvolutionarySolution<FT>(sol1);
		// Find out range to crossover
		//cout << "New call of crossingOver";

		//std::random_device seeder_mod;
		//std::mt19937 engine(seeder_mod());

		//std::uniform_int_distribution<int> dist_dim_start(0, (sol1._dim-2)/2);
		//int start_dim = dist_dim_start(engine);
		int start_dim = (rand()%(sol1._dim-2)/2);

		//std::uniform_int_distribution<int> dist_dim_end(start_dim+1, (sol1._dim-1)/2);
		//int end_dim  = dist_dim_end(engine);
		int end_dim  = std::max(rand () % ((sol1._dim-1)/2 - start_dim) + start_dim + 1, min(start_dim+8, (sol1._dim-2)/2));
		//cout << start_dim << " // " << end_dim << endl;



		if(start_dim >= end_dim || end_dim*2+1 >= sol1._dim) {
			int stop=0;
			cerr << "EvolutionarySolutionCrossing over with start " << start_dim
					<< " and end " << end_dim << endl;
			cin >> stop;
		}

		// Copy the values
		for(int i = start_dim*2; i <= end_dim*2+1; i++) {
			sol_res._prun_func[i] = sol2._prun_func[i];
		}

		// Invalid solutions have to be repaired
		FT low_crossval = sol_res._prun_func[start_dim*2];
		FT high_crossval = sol_res._prun_func[end_dim*2];

		// Check below
		for(int i=start_dim*2; i>=1; i--) {
			if(low_crossval < sol_res._prun_func[i-1]) {
				sol_res._prun_func[i-1] = low_crossval;
			}
			else {
				break;
			}
		}

		// Check above
		for(int i=end_dim*2+1; i <= sol_res._dim-3; i++) {
			if(high_crossval > sol_res._prun_func[i+1]) {
				sol_res._prun_func[i+1] = high_crossval;
			}
			else {
				break;
			}
		}

		// Check feasiability
		/*for(long unsigned int i=0; i<sol_res._prun_func.size()-1; i++) {
				if(sol_res._prun_func[i] > sol_res._prun_func[i+1]) {
					cout << "Prunfunc invalid! :" << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
					<< i << " (" << sol_res._prun_func[i] << " / " << sol_res._prun_func[i+1] << ")" << endl << std::flush;


					cout << "EvolutionarySolutionCrossing over with start " << start_dim*2
							<< " and end " << end_dim*2
							<< " with ID(" << sol1.getRandomID() << ") and ID("
							<< sol2.getRandomID() << ")"
							<< endl << std::flush;

					// Copy the values
					cout << "Start writing: " << endl;
					for(int i = start_dim*2; i <= end_dim*2+1; i++) {
						sol_res._prun_func[i] = sol2._prun_func[i];
						cout << "[" << i << "]: " << sol2._prun_func[i] << endl;
					}
					cout << "End writing" << endl;

					sol_res.printSolution(true);

					// Check below
					cout << "Below: ";
					for(int i=start_dim*2; i>=1; i--) {
						cout << low_crossval << "/" << sol_res._prun_func[i-1] << " ";
						if(low_crossval < sol_res._prun_func[i-1]) {
							cout << i << endl;
						}
						else {
							break;
						}
					}

					cout << endl << "Above: ";
					// Check above
					for(int i=end_dim*2+1; i <= sol_res._dim-3; i++) {
						cout << high_crossval << "/" << sol_res._prun_func[i+1] << " ";
						if(high_crossval > sol_res._prun_func[i+1]) {
							cout << i << endl;
						}
						else {
							break;
						}
					}
					cout << endl;


					exit(-1000);
				}
		}*/
		//cout << "Test passed" << endl;

		// Recalc probs and costs!
		sol_res.resetProbs();
		sol_res.calculateCosts();
		return sol_res;
	}

};




#endif /* SRC_EVOLUTIONPOPULATION_HPP_ */
