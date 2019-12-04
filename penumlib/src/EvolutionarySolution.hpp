/*
 * EvolutionarySolution.hpp
 *
 *  Created on: 15.10.2019
 *      Author: mburger
 */

#ifndef SRC_EVOLUTIONARYSOLUTION_HPP_
#define SRC_EVOLUTIONARYSOLUTION_HPP_

#include "AnnealingSolution.hpp"

template <class FT1>
class EvolutionarySolution : public AnnealingSolution<FT1> {
public:
	EvolutionarySolution() : AnnealingSolution<FT1>(){

	}

	EvolutionarySolution(BKZBenchmarker<FT1>* bench_mark, const vector<FT1>& prun_func,
			const AnnealInfo<FT1>& ainfo) : AnnealingSolution<FT1> (bench_mark, prun_func, ainfo) {

	}

	EvolutionarySolution(EvolutionarySolution& other) : AnnealingSolution<FT1>(other) { //_bench_mark (other._bench_mark)
		//this->_bench_mark = other._bench_mark;
		/*this->_costs = other._costs;
		this->_ainfo = other._ainfo;
		this->_dim = other._dim;

		this->_succ_prob.reserve(other._ainfo._number_of_random_bases);
		this->_t_nodes.reserve(other._ainfo._number_of_random_bases);
		this->_est_nodes.reserve(other._ainfo._number_of_random_bases);
		this->_t_reduction.reserve(other._ainfo._number_of_random_bases);
		this->_t_enums.reserve(other._ainfo._number_of_random_bases);
		this->_base_costs.reserve(other._ainfo._number_of_random_bases);

		this->_succ_prob.resize(other._ainfo._number_of_random_bases);
		this->_t_nodes.resize(other._ainfo._number_of_random_bases);
		this->_est_nodes.resize(other._ainfo._number_of_random_bases);
		this->_t_reduction.resize(other._ainfo._number_of_random_bases);
		this->_t_enums.resize(other._ainfo._number_of_random_bases);
		this->_base_costs.resize(other._ainfo._number_of_random_bases);

		this->_probs_initialized.reserve(other._ainfo._number_of_random_bases);
		this->_probs_initialized.resize(other._ainfo._number_of_random_bases);
		this->_probs.reserve(other._ainfo._number_of_random_bases);
		this->_probs.resize(other._ainfo._number_of_random_bases);

		for(int i=0; i < other._ainfo._number_of_random_bases; i++) {
			this->_probs_initialized[i] = other._probs_initialized[i];
			this->_succ_prob[i] = other._succ_prob[i];
			this->_t_nodes[i] = other._t_nodes[i];
			this->_est_nodes[i] = other._est_nodes[i];
			this->_t_reduction[i] = other._t_reduction[i];
			this->_t_enums[i] = other._t_enums[i];
			this->_base_costs[i] = other._base_costs[i];
			this->_probs[i].clear();

			this->_probs[i].reserve(this->_dim);
			for(auto it = other._probs[i].begin(); it != other._probs[i].end(); it++) {
				this->_probs[i].push_back(*it);
			}
		}

		this->_costs_calculated = other._costs_calculated;

		this->_prun_func.clear();
		for(auto it = other._prun_func.begin(); it != other._prun_func.end(); it++) {
			this->_prun_func.push_back(*it);
		}

		this->_bkz_info.betas.first = other._bkz_info.betas.first;
		this->_bkz_info.betas.second = other._bkz_info.betas.second;
		this->_bkz_info.rtime.reserve(other._ainfo._number_of_random_bases);
		this->_bkz_info.rtime.resize(other._ainfo._number_of_random_bases);
		this->_bkz_info.seed.reserve(other._ainfo._number_of_random_bases);
		this->_bkz_info.seed.resize(other._ainfo._number_of_random_bases);
		this->_bkz_info.amax.reserve(other._ainfo._number_of_random_bases);
		this->_bkz_info.amax.resize(other._ainfo._number_of_random_bases);

		for(int i=0; i < this->_ainfo._number_of_random_bases; i++) {
			this->_bkz_info.rtime[i] = other._bkz_info.rtime[i];
			this->_bkz_info.seed[i] = other._bkz_info.seed[i];
			this->_bkz_info.amax[i] = other._bkz_info.amax[i];
		}

		this->setRandomID();*/
	}

	EvolutionarySolution(const EvolutionarySolution& other) : AnnealingSolution<FT1> (other) {
		//this->_bench_mark = other._bench_mark;
		/*this->_costs = other._costs;
		this->_ainfo = other._ainfo;
		this->_dim = other._dim;

		this->_succ_prob.reserve(other._ainfo._number_of_random_bases);
		this->_t_nodes.reserve(other._ainfo._number_of_random_bases);
		this->_est_nodes.reserve(other._ainfo._number_of_random_bases);
		this->_t_reduction.reserve(other._ainfo._number_of_random_bases);
		this->_t_enums.reserve(other._ainfo._number_of_random_bases);
		this->_base_costs.reserve(other._ainfo._number_of_random_bases);

		this->_succ_prob.resize(other._ainfo._number_of_random_bases);
		this->_t_nodes.resize(other._ainfo._number_of_random_bases);
		this->_est_nodes.resize(other._ainfo._number_of_random_bases);
		this->_t_reduction.resize(other._ainfo._number_of_random_bases);
		this->_t_enums.resize(other._ainfo._number_of_random_bases);
		this->_base_costs.resize(other._ainfo._number_of_random_bases);

		this->_probs_initialized.reserve(other._ainfo._number_of_random_bases);
		this->_probs_initialized.resize(other._ainfo._number_of_random_bases);
		this->_probs.reserve(other._ainfo._number_of_random_bases);
		this->_probs.resize(other._ainfo._number_of_random_bases);

		for(int i=0; i < other._ainfo._number_of_random_bases; i++) {
			this->_probs_initialized[i] = other._probs_initialized[i];
			this->_succ_prob[i] = other._succ_prob[i];
			this->_t_nodes[i] = other._t_nodes[i];
			this->_est_nodes[i] = other._est_nodes[i];
			this->_t_reduction[i] = other._t_reduction[i];
			this->_t_enums[i] = other._t_enums[i];
			this->_base_costs[i] = other._base_costs[i];
			this->_probs[i].clear();

			this->_probs[i].reserve(this->_dim);
			for(auto it = other._probs[i].begin(); it != other._probs[i].end(); it++) {
				this->_probs[i].push_back(*it);
			}
		}

		this->_costs_calculated = other._costs_calculated;

		this->_prun_func.clear();
		for(auto it = other._prun_func.begin(); it != other._prun_func.end(); it++) {
			this->_prun_func.push_back(*it);
		}

		this->_bkz_info.betas.first = other._bkz_info.betas.first;
		this->_bkz_info.betas.second = other._bkz_info.betas.second;
		this->_bkz_info.rtime.reserve(other._ainfo._number_of_random_bases);
		this->_bkz_info.rtime.resize(other._ainfo._number_of_random_bases);
		this->_bkz_info.seed.reserve(other._ainfo._number_of_random_bases);
		this->_bkz_info.seed.resize(other._ainfo._number_of_random_bases);
		this->_bkz_info.amax.reserve(other._ainfo._number_of_random_bases);
		this->_bkz_info.amax.resize(other._ainfo._number_of_random_bases);

		for(int i=0; i < this->_ainfo._number_of_random_bases; i++) {
			this->_bkz_info.rtime[i] = other._bkz_info.rtime[i];
			this->_bkz_info.seed[i] = other._bkz_info.seed[i];
			this->_bkz_info.amax[i] = other._bkz_info.amax[i];
		}

		this->setRandomID();*/
	}

	EvolutionarySolution& operator= (const EvolutionarySolution & other) {
		//this->_bench_mark = other._bench_mark;
		this->_costs = other._costs;
		this->_ainfo = other._ainfo;
		this->_dim = other._dim;

		this->_succ_prob.reserve(other._ainfo._number_of_random_bases);
		this->_prob_thread_shots.reserve(other._ainfo._number_of_random_bases);
		this->_no_of_shots.reserve(other._ainfo._number_of_random_bases);
		this->_t_nodes.reserve(other._ainfo._number_of_random_bases);
		this->_est_nodes.reserve(other._ainfo._number_of_random_bases);
		this->_t_reduction.reserve(other._ainfo._number_of_random_bases);
		this->_t_enums.reserve(other._ainfo._number_of_random_bases);
		this->_base_costs.reserve(other._ainfo._number_of_random_bases);

		this->_succ_prob.resize(other._ainfo._number_of_random_bases);
		this->_prob_thread_shots.resize(other._ainfo._number_of_random_bases);
		this->_no_of_shots.resize(other._ainfo._number_of_random_bases);
		this->_t_nodes.resize(other._ainfo._number_of_random_bases);
		this->_est_nodes.resize(other._ainfo._number_of_random_bases);
		this->_t_reduction.resize(other._ainfo._number_of_random_bases);
		this->_t_enums.resize(other._ainfo._number_of_random_bases);
		this->_base_costs.resize(other._ainfo._number_of_random_bases);

		this->_probs_initialized.reserve(other._ainfo._number_of_random_bases);
		this->_probs_initialized.resize(other._ainfo._number_of_random_bases);
		this->_probs.reserve(other._ainfo._number_of_random_bases);
		this->_probs.resize(other._ainfo._number_of_random_bases);

		for(int i=0; i < other._ainfo._number_of_random_bases; i++) {
			this->_probs_initialized[i] = other._probs_initialized[i];
			this->_prob_thread_shots[i] = other._prob_thread_shots[i];
			this->_no_of_shots[i] = other._no_of_shots[i];
			this->_succ_prob[i] = other._succ_prob[i];
			this->_t_nodes[i] = other._t_nodes[i];
			this->_est_nodes[i] = other._est_nodes[i];
			this->_t_reduction[i] = other._t_reduction[i];
			this->_t_enums[i] = other._t_enums[i];
			this->_base_costs[i] = other._base_costs[i];
			this->_probs[i].clear();
			for(auto it = other._probs[i].begin(); it != other._probs[i].end(); it++) {
				this->_probs[i].push_back(*it);
			}
		}

		this->_costs_calculated = other._costs_calculated;

		this->_prun_func.clear();
		for(auto it = other._prun_func.begin(); it != other._prun_func.end(); it++) {
			this->_prun_func.push_back(*it);
		}

		this->_bkz_info.betas.first = other._bkz_info.betas.first;
		this->_bkz_info.betas.second = other._bkz_info.betas.second ;
		this->_bkz_info.rtime.reserve(other._ainfo._number_of_random_bases);
		this->_bkz_info.rtime.resize(other._ainfo._number_of_random_bases);
		this->_bkz_info.seed.reserve(other._ainfo._number_of_random_bases);
		this->_bkz_info.seed.resize(other._ainfo._number_of_random_bases);
		this->_bkz_info.amax.reserve(other._ainfo._number_of_random_bases);
		this->_bkz_info.amax.resize(other._ainfo._number_of_random_bases);

		for(int i=0; i < this->_ainfo._number_of_random_bases; i++) {
			this->_bkz_info.rtime[i] = other._bkz_info.rtime[i];
			this->_bkz_info.seed[i] = other._bkz_info.seed[i];
			this->_bkz_info.amax[i] = other._bkz_info.amax[i];
		}

		//this->_id = other._id;
		this->setRandomID();
		return *this;
	}

	EvolutionarySolution& operator= (EvolutionarySolution & other) {
			//this->_bench_mark = other._bench_mark;
			this->_costs = other._costs;
			this->_ainfo = other._ainfo;
			this->_dim = other._dim;

			this->_succ_prob.reserve(other._ainfo._number_of_random_bases);
			this->_prob_thread_shots.reserve(other._ainfo._number_of_random_bases);
			this->_no_of_shots.reserve(other._ainfo._number_of_random_bases);
			this->_t_nodes.reserve(other._ainfo._number_of_random_bases);
			this->_est_nodes.reserve(other._ainfo._number_of_random_bases);
			this->_t_reduction.reserve(other._ainfo._number_of_random_bases);
			this->_t_enums.reserve(other._ainfo._number_of_random_bases);
			this->_base_costs.reserve(other._ainfo._number_of_random_bases);

			this->_succ_prob.resize(other._ainfo._number_of_random_bases);
			this->_prob_thread_shots.resize(other._ainfo._number_of_random_bases);
			this->_no_of_shots.resize(other._ainfo._number_of_random_bases);
			this->_t_nodes.resize(other._ainfo._number_of_random_bases);
			this->_est_nodes.resize(other._ainfo._number_of_random_bases);
			this->_t_reduction.resize(other._ainfo._number_of_random_bases);
			this->_t_enums.resize(other._ainfo._number_of_random_bases);
			this->_base_costs.resize(other._ainfo._number_of_random_bases);

			this->_probs_initialized.reserve(other._ainfo._number_of_random_bases);
			this->_probs_initialized.resize(other._ainfo._number_of_random_bases);
			this->_probs.reserve(other._ainfo._number_of_random_bases);
			this->_probs.resize(other._ainfo._number_of_random_bases);

			for(int i=0; i < other._ainfo._number_of_random_bases; i++) {
				this->_probs_initialized[i] = other._probs_initialized[i];
				this->_succ_prob[i] = other._succ_prob[i];
				this->_prob_thread_shots[i] = other._prob_thread_shots[i];
				this->_no_of_shots[i] = other._no_of_shots[i];
				this->_t_nodes[i] = other._t_nodes[i];
				this->_est_nodes[i] = other._est_nodes[i];
				this->_t_reduction[i] = other._t_reduction[i];
				this->_t_enums[i] = other._t_enums[i];
				this->_base_costs[i] = other._base_costs[i];
				this->_probs[i].clear();
				for(auto it = other._probs[i].begin(); it != other._probs[i].end(); it++) {
					this->_probs[i].push_back(*it);
				}
			}

			this->_costs_calculated = other._costs_calculated;

			this->_prun_func.clear();
			for(auto it = other._prun_func.begin(); it != other._prun_func.end(); it++) {
				this->_prun_func.push_back(*it);
			}

			this->_bkz_info.betas.first = other._bkz_info.betas.first;
			this->_bkz_info.betas.second = other._bkz_info.betas.second ;
			this->_bkz_info.rtime.reserve(other._ainfo._number_of_random_bases);
			this->_bkz_info.rtime.resize(other._ainfo._number_of_random_bases);
			this->_bkz_info.seed.reserve(other._ainfo._number_of_random_bases);
			this->_bkz_info.seed.resize(other._ainfo._number_of_random_bases);
			this->_bkz_info.amax.reserve(other._ainfo._number_of_random_bases);
			this->_bkz_info.amax.resize(other._ainfo._number_of_random_bases);

			for(int i=0; i < this->_ainfo._number_of_random_bases; i++) {
				this->_bkz_info.rtime[i] = other._bkz_info.rtime[i];
				this->_bkz_info.seed[i] = other._bkz_info.seed[i];
				this->_bkz_info.amax[i] = other._bkz_info.amax[i];
			}

			//this->_id = other._id;
			this->setRandomID();
			//cout << "Changing ID from " << other._id << " to " << this->_id << endl;
			return *this;
		}



	/** Checks whether the solution is valid **/
	bool checkValidity() {
		int dim_half = this->_dim / 2;

		for (int i=0; i < dim_half; i++) {
			// Neighbors even-uneven have to be equal
			if(this->_prun_func[i*2] != this->_prun_func[i*2 + 1])
				return false;

			// Previous step must not be bigger
			if(i > 0 && this->_prun_func[(i-1)*2] > this->_prun_func[i*2])
				return false;

			// Next step must be higher
			if(i <  (dim_half -1) && this->_prun_func[i*2] > this->_prun_func[(i+1)*2])
				return false;

			if(this->_prun_func[i*2] > 1)
				return false;
		}
		return true;
	}

	int findMaxRecombineEntry(const EvolutionarySolution& sol1, const EvolutionarySolution& sol2) {
		int i=0;
		for(; i < sol1._dim/2; i++) {
			if(sol1._prun_func[i] > sol2._prun_func[i])
				break;
		}
		return i-1;
	}

	FT1 getFitness() const {
		return (FT1) ((1000000.0 / this->_costs) * (1000000.0 / this->_costs));
	}

	//void setFitness(FT1 f) {
	//	this->_fitness = f;
	//}

	template <class FT>
	friend	EvolutionarySolution<FT> recombine(const EvolutionarySolution<FT>& sol1, const EvolutionarySolution<FT>& sol2);
	template <class FT>
	friend	EvolutionarySolution<FT> crossingOver(const EvolutionarySolution<FT>& sol1, const EvolutionarySolution<FT>& sol2);
	template <class FT>
	friend	EvolutionarySolution<FT> shift(const EvolutionarySolution<FT>& sol1);

private:
	//BKZBenchmarker<FT1>* const _bench_mark;
	//FT1 _fitness;
};



#endif /* SRC_EVOLUTIONARYSOLUTION_HPP_ */
