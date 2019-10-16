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


using namespace std;
//template <class FT>
//using individual = std::pair <int, EvolutionarySolution<FT>> ;

template <class FT>
struct classcomp {
	bool operator() (const EvolutionarySolution<FT>& lhs, const EvolutionarySolution<FT>& rhs) const
	{
		return lhs.getFitness() > rhs.getFitness();
		//return lhs.getRandomID() > rhs.getRandomID();
	}
};

template <class FT>
class EvolutionPopulation {
public:
	EvolutionPopulation(int target_size) {
		this->_population.clear();
		this->_target_size = target_size;
		cout << "Setting target size to " << this->_target_size << endl;
		this->_best_individual = EvolutionarySolution<FT>();
		this->_max_costs = 0;
	}

	EvolutionPopulation() {
		this->_population.clear();
		this->_target_size = 101;
		cout << "Setting default target size to " << this->_target_size << endl;
		this->_best_individual = EvolutionarySolution<FT>();
		this->_max_costs = 0;
	}

	~EvolutionPopulation () {

	}

	EvolutionPopulation(EvolutionPopulation& other) {
		this->_target_size = other._target_size;
		this->_max_costs = other._max_costs;
		this->_best_individual = other._best_individual;
		this->_population = other._population;
	}

	EvolutionPopulation(const EvolutionPopulation& other) {
		this->_target_size = other._target_size;
		this->_max_costs = other._max_costs;
		this->_best_individual = other._best_individual;
		this->_population = other._population;
	}

	EvolutionPopulation& operator= (const EvolutionPopulation & other) {
		_target_size = other._target_size;
		_max_costs = other._max_costs;
		_best_individual = other._best_individual;
		_population = other._population;
		return *this;
	}

	int insertIndividual (EvolutionarySolution<FT>& indi) {
		EvolutionarySolution<FT> indi_copy = indi;
		this->_population.insert(indi_copy);
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

		for(int i=0; i<_target_size/1; i++) {
			//cout << i << "/" << _target_size << endl;
			parents.insert( selectionWithRoullette() );
		}

		_population = parents;


		return 0;
	}

	void printCostsOfIndividuals() {
		auto iti = _population.begin();

		for(;iti != _population.end(); ++iti) {
			cout << "ID(" << iti->getRandomID() << ")"
			<< iti->getFitness() << endl;
		}

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

	EvolutionarySolution<FT> selectionWithRoullette() {
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

protected:
	std::set<EvolutionarySolution<FT>, classcomp<FT>> _population; /**/
	EvolutionarySolution<FT> _best_individual;
	double _max_costs;
	int _target_size;

};




#endif /* SRC_EVOLUTIONPOPULATION_HPP_ */
