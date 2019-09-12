/*
 * AbstractAnnealingSolution.hpp
 *
 *  Created on: 26.02.2019
 *      Author: mburger
 */

#ifndef SRC_ANNEALINGSOLUTION_HPP_
#define SRC_ANNEALINGSOLUTION_HPP_

#include <random>
#include <vector>
#include "BKZBenchmarker.hpp"
#include "extremePruningFunction.hpp"
#include <limits>

using namespace std;

template <class FT1>
class AnnealingSolution {
public:
	AnnealingSolution()
	{
		setRandomID();
		_costs_calculated = false;
		_ainfo = AnnealInfo<FT1>();
	}

	AnnealingSolution(BKZBenchmarker<FT1>* bench_mark, const vector<FT1>& prun_func, const AnnealInfo<FT1>& ainfo) {
		setRandomID();

		_costs_calculated = false;
		_ainfo = ainfo;
		_dim = prun_func.size();

		_succ_prob.reserve(ainfo._number_of_random_bases);
		_t_nodes.reserve(ainfo._number_of_random_bases);
		_est_nodes.reserve(ainfo._number_of_random_bases);
		_t_reduction.reserve(ainfo._number_of_random_bases);
		_t_enums.reserve(ainfo._number_of_random_bases);
		_base_costs.reserve(ainfo._number_of_random_bases);
		_probs_initialized.reserve(ainfo._number_of_random_bases);

		_succ_prob.resize(ainfo._number_of_random_bases);
		_t_nodes.resize(ainfo._number_of_random_bases);
		_est_nodes.resize(ainfo._number_of_random_bases);
		_t_reduction.resize(ainfo._number_of_random_bases);
		_t_enums.resize(ainfo._number_of_random_bases);
		_base_costs.resize(ainfo._number_of_random_bases);
		_probs_initialized.resize(ainfo._number_of_random_bases);

		_probs.reserve(ainfo._number_of_random_bases);
		_probs.resize(ainfo._number_of_random_bases);

		for(int i=0; i < ainfo._number_of_random_bases; i++) {
			_probs[i].reserve(_dim);
			_probs[i].resize(_dim);
		}

		_bench_mark = bench_mark;
		for(auto it=prun_func.begin(); it!=prun_func.end();it++) {
			_prun_func.push_back(*it);
		}

		BKZInfo<FT1>& bkz_info = bench_mark->next();

		_bkz_info.betas.first = bkz_info.betas.first;
		_bkz_info.betas.second = bkz_info.betas.second;
		_bkz_info.rtime.reserve(_ainfo._number_of_random_bases);
		_bkz_info.rtime.resize(_ainfo._number_of_random_bases);
		_bkz_info.seed.reserve(_ainfo._number_of_random_bases);
		_bkz_info.seed.resize(_ainfo._number_of_random_bases);
		_bkz_info.amax.reserve(_ainfo._number_of_random_bases);
		_bkz_info.amax.resize(_ainfo._number_of_random_bases);

		for(int i=0; i < _ainfo._number_of_random_bases; i++) {
			_probs_initialized[i] = false;
			_bkz_info.rtime[i] = bkz_info.rtime[i];
			_bkz_info.seed[i] = bkz_info.seed[i];
			_bkz_info.amax[i] = bkz_info.amax[i];
		}
		_costs_calculated = true;

		setOpimtalBetaPair();

		for(int i=0; i < _ainfo._number_of_random_bases; i++) {
			_probs_initialized[i] = true;
		}
	}

	AnnealingSolution(AnnealingSolution& other) {
		_bench_mark = other._bench_mark;
		_costs = other._costs;
		_ainfo = other._ainfo;
		_dim = other._dim;

		_succ_prob.reserve(other._ainfo._number_of_random_bases);
		_t_nodes.reserve(other._ainfo._number_of_random_bases);
		_est_nodes.reserve(other._ainfo._number_of_random_bases);
		_t_reduction.reserve(other._ainfo._number_of_random_bases);
		_t_enums.reserve(other._ainfo._number_of_random_bases);
		_base_costs.reserve(other._ainfo._number_of_random_bases);

		_succ_prob.resize(other._ainfo._number_of_random_bases);
		_t_nodes.resize(other._ainfo._number_of_random_bases);
		_est_nodes.resize(other._ainfo._number_of_random_bases);
		_t_reduction.resize(other._ainfo._number_of_random_bases);
		_t_enums.resize(other._ainfo._number_of_random_bases);
		_base_costs.resize(other._ainfo._number_of_random_bases);

		_probs_initialized.reserve(other._ainfo._number_of_random_bases);
		_probs_initialized.resize(other._ainfo._number_of_random_bases);
		_probs.reserve(other._ainfo._number_of_random_bases);
		_probs.resize(other._ainfo._number_of_random_bases);

		for(int i=0; i < other._ainfo._number_of_random_bases; i++) {
			_probs_initialized[i] = other._probs_initialized[i];
			_succ_prob[i] = other._succ_prob[i];
			_t_nodes[i] = other._t_nodes[i];
			_est_nodes[i] = other._est_nodes[i];
			_t_reduction[i] = other._t_reduction[i];
			_t_enums[i] = other._t_enums[i];
			_base_costs[i] = other._base_costs[i];
			this->_probs[i].clear();

			_probs[i].reserve(_dim);
			for(auto it = other._probs[i].begin(); it != other._probs[i].end(); it++) {
				this->_probs[i].push_back(*it);
			}
		}

		_costs_calculated = other._costs_calculated;

		this->_prun_func.clear();
		for(auto it = other._prun_func.begin(); it != other._prun_func.end(); it++) {
			this->_prun_func.push_back(*it);
		}

		_bkz_info.betas.first = other._bkz_info.betas.first;
		_bkz_info.betas.second = other._bkz_info.betas.second;
		_bkz_info.rtime.reserve(other._ainfo._number_of_random_bases);
		_bkz_info.rtime.resize(other._ainfo._number_of_random_bases);
		_bkz_info.seed.reserve(other._ainfo._number_of_random_bases);
		_bkz_info.seed.resize(other._ainfo._number_of_random_bases);
		_bkz_info.amax.reserve(other._ainfo._number_of_random_bases);
		_bkz_info.amax.resize(other._ainfo._number_of_random_bases);

		for(int i=0; i < _ainfo._number_of_random_bases; i++) {
			_bkz_info.rtime[i] = other._bkz_info.rtime[i];
			_bkz_info.seed[i] = other._bkz_info.seed[i];
			_bkz_info.amax[i] = other._bkz_info.amax[i];
		}

		setRandomID();
	}

	AnnealingSolution(const AnnealingSolution& other) {
		_bench_mark = other._bench_mark;
		_costs = other._costs;
		_ainfo = other._ainfo;
		_dim = other._dim;

		_succ_prob.reserve(other._ainfo._number_of_random_bases);
		_t_nodes.reserve(other._ainfo._number_of_random_bases);
		_est_nodes.reserve(other._ainfo._number_of_random_bases);
		_t_reduction.reserve(other._ainfo._number_of_random_bases);
		_t_enums.reserve(other._ainfo._number_of_random_bases);
		_base_costs.reserve(other._ainfo._number_of_random_bases);

		_succ_prob.resize(other._ainfo._number_of_random_bases);
		_t_nodes.resize(other._ainfo._number_of_random_bases);
		_est_nodes.resize(other._ainfo._number_of_random_bases);
		_t_reduction.resize(other._ainfo._number_of_random_bases);
		_t_enums.resize(other._ainfo._number_of_random_bases);
		_base_costs.resize(other._ainfo._number_of_random_bases);

		_probs_initialized.reserve(other._ainfo._number_of_random_bases);
		_probs_initialized.resize(other._ainfo._number_of_random_bases);
		_probs.reserve(other._ainfo._number_of_random_bases);
		_probs.resize(other._ainfo._number_of_random_bases);

		for(int i=0; i < other._ainfo._number_of_random_bases; i++) {
			_probs_initialized[i] = other._probs_initialized[i];
			_succ_prob[i] = other._succ_prob[i];
			_t_nodes[i] = other._t_nodes[i];
			_est_nodes[i] = other._est_nodes[i];
			_t_reduction[i] = other._t_reduction[i];
			_t_enums[i] = other._t_enums[i];
			_base_costs[i] = other._base_costs[i];
			this->_probs[i].clear();

			_probs[i].reserve(_dim);
			for(auto it = other._probs[i].begin(); it != other._probs[i].end(); it++) {
				this->_probs[i].push_back(*it);
			}
		}

		_costs_calculated = other._costs_calculated;

		this->_prun_func.clear();
		for(auto it = other._prun_func.begin(); it != other._prun_func.end(); it++) {
			this->_prun_func.push_back(*it);
		}

		_bkz_info.betas.first = other._bkz_info.betas.first;
		_bkz_info.betas.second = other._bkz_info.betas.second;
		_bkz_info.rtime.reserve(other._ainfo._number_of_random_bases);
		_bkz_info.rtime.resize(other._ainfo._number_of_random_bases);
		_bkz_info.seed.reserve(other._ainfo._number_of_random_bases);
		_bkz_info.seed.resize(other._ainfo._number_of_random_bases);
		_bkz_info.amax.reserve(other._ainfo._number_of_random_bases);
		_bkz_info.amax.resize(other._ainfo._number_of_random_bases);

		for(int i=0; i < _ainfo._number_of_random_bases; i++) {
			_bkz_info.rtime[i] = other._bkz_info.rtime[i];
			_bkz_info.seed[i] = other._bkz_info.seed[i];
			_bkz_info.amax[i] = other._bkz_info.amax[i];
		}

		setRandomID();
	}

	AnnealingSolution& operator= (const AnnealingSolution & other) {
		_bench_mark = other._bench_mark;
		_costs = other._costs;
		_ainfo = other._ainfo;
		_dim = other._dim;

		_succ_prob.reserve(other._ainfo._number_of_random_bases);
		_t_nodes.reserve(other._ainfo._number_of_random_bases);
		_est_nodes.reserve(other._ainfo._number_of_random_bases);
		_t_reduction.reserve(other._ainfo._number_of_random_bases);
		_t_enums.reserve(other._ainfo._number_of_random_bases);
		_base_costs.reserve(other._ainfo._number_of_random_bases);

		_succ_prob.resize(other._ainfo._number_of_random_bases);
		_t_nodes.resize(other._ainfo._number_of_random_bases);
		_est_nodes.resize(other._ainfo._number_of_random_bases);
		_t_reduction.resize(other._ainfo._number_of_random_bases);
		_t_enums.resize(other._ainfo._number_of_random_bases);
		_base_costs.resize(other._ainfo._number_of_random_bases);

		_probs_initialized.reserve(other._ainfo._number_of_random_bases);
		_probs_initialized.resize(other._ainfo._number_of_random_bases);
		_probs.reserve(other._ainfo._number_of_random_bases);
		_probs.resize(other._ainfo._number_of_random_bases);

		for(int i=0; i < other._ainfo._number_of_random_bases; i++) {
			_probs_initialized[i] = other._probs_initialized[i];
			_succ_prob[i] = other._succ_prob[i];
			_t_nodes[i] = other._t_nodes[i];
			_est_nodes[i] = other._est_nodes[i];
			_t_reduction[i] = other._t_reduction[i];
			_t_enums[i] = other._t_enums[i];
			_base_costs[i] = other._base_costs[i];
			this->_probs[i].clear();
			for(auto it = other._probs[i].begin(); it != other._probs[i].end(); it++) {
				this->_probs[i].push_back(*it);
			}
		}

		_costs_calculated = other._costs_calculated;

		this->_prun_func.clear();
		for(auto it = other._prun_func.begin(); it != other._prun_func.end(); it++) {
			this->_prun_func.push_back(*it);
		}

		_bkz_info.betas.first = other._bkz_info.betas.first;
		_bkz_info.betas.second = other._bkz_info.betas.second ;
		_bkz_info.rtime.reserve(other._ainfo._number_of_random_bases);
		_bkz_info.rtime.resize(other._ainfo._number_of_random_bases);
		_bkz_info.seed.reserve(other._ainfo._number_of_random_bases);
		_bkz_info.seed.resize(other._ainfo._number_of_random_bases);
		_bkz_info.amax.reserve(other._ainfo._number_of_random_bases);
		_bkz_info.amax.resize(other._ainfo._number_of_random_bases);

		for(int i=0; i < _ainfo._number_of_random_bases; i++) {
			_bkz_info.rtime[i] = other._bkz_info.rtime[i];
			_bkz_info.seed[i] = other._bkz_info.seed[i];
			_bkz_info.amax[i] = other._bkz_info.amax[i];
		}

		setRandomID();
		return *this;
	}

	~AnnealingSolution() {
	}

	void printSolutionToSStream(stringstream& ss)
	{
		ss << "PreBeta: " << this->_bkz_info.betas.first << " / " << "Beta: " << this->_bkz_info.betas.second << endl;

		ss << "Time BKZ: [ ";
		for(int i=0; i < _ainfo._number_of_random_bases; i++) {
			ss << this->_t_reduction[i] << "s ";
		}
		ss << "]" << endl;

		ss << "#Nodes: [ ";
		for(int i=0; i < _ainfo._number_of_random_bases; i++) {
			ss << this->_est_nodes[i] << " ";
		}
		ss << "]" << endl;

		ss << "Time ENUM: [ ";
		for(int i=0; i < _ainfo._number_of_random_bases; i++) {
			ss << this->_t_enums[i] << "s ";
		}
		ss << "]" << endl;

		ss << "Probs: [ ";
		for(int i=0; i < _ainfo._number_of_random_bases; i++) {
			ss << this->_succ_prob[i] * 100 << "% ";
		}
		ss << "]" << endl;

		ss << "Total costs: [ ";
		for(int i=0; i < _ainfo._number_of_random_bases; i++) {
			ss << this->_base_costs[i] << "s ";
		}
		ss << "]" << endl;

		// Print the probability values
		for(auto it=_prun_func.rbegin(); it != _prun_func.rend(); ++it) {
			ss << *it << endl;
		}
	}

	void printSolution() {
		cout << "PreBeta: " << this->_bkz_info.betas.first << " / " << "Beta: " << this->_bkz_info.betas.second << endl;

		cout << "Time BKZ: [ ";
		for(int i=0; i < _ainfo._number_of_random_bases; i++) {
			cout << this->_t_reduction[i] << "s ";
		}
		cout << "]" << endl;

		cout << "#Nodes: [ ";
		for(int i=0; i < _ainfo._number_of_random_bases; i++) {
			cout << this->_est_nodes[i] << " ";
		}
		cout << "]" << endl;

		cout << "Time ENUM: [ ";
		for(int i=0; i < _ainfo._number_of_random_bases; i++) {
			cout << this->_t_enums[i] << "s ";
		}
		cout << "]" << endl;

		cout << "Probs: [ ";
		for(int i=0; i < _ainfo._number_of_random_bases; i++) {
			cout << this->_succ_prob[i] * 100 << "% ";
		}
		cout << "]" << endl;

		cout << "Total costs: [ ";
		for(int i=0; i < _ainfo._number_of_random_bases; i++) {
			cout << this->_base_costs[i] << "s ";
		}
		cout << "]" << endl;

		// Print the probability values
		for(auto it=_prun_func.rbegin(); it != _prun_func.rend(); ++it) {
			cout << *it << endl;
		}

		cout << "Overall cost: " << _costs << endl;
	}

	FT1 getCost() {
		return _costs;
	}

	int setToNumberedBetaPair (unsigned int no) {
		std::map<BetaPair, BKZInfo<FT1>>& the_map = this->_bench_mark->getDataMap();

		if(no < 0) {
			cerr << "Requested negative than BKZ-Map size. Setting to best pair.";
			setOpimtalBetaPair();
			return -1;
		}
        unsigned int length = the_map.size();
		unsigned int ass_num = no % length;

		//
        
		if(no >= length) {
			unsigned int middle = (int)floor(the_map.size()/2.0);
            
            if(length%2 != 0) {
                if(ass_num == middle) {
                    // Do nothing because mapping already correct
                }
                
                else if(ass_num%2 != 0) {
                    ass_num = middle - (int)(ceil(ass_num / 2.0));
                }
                
                else  {
                    ass_num = middle + (int)(ceil(ass_num / 2.0));
                }
            }
            
            // even length
            else {
                cerr << "Assigment of threads to even number of beta-configs not implemented yet!" << endl;
                ass_num = no % length;
            }
		}

		// Find the numbered entry
		auto it = the_map.begin();
		for(unsigned int i=0; i<ass_num;i++) {
			++it;
		}

		this->_bkz_info = it->second;
        //#pragma omp critical
        //cout << "Assiging " << _bkz_info.betas.second << " to " << no << " with ass_num " << ass_num  << endl << std::flush;

		// If beta changed, the costs also changed
		this->calculateCosts();
		return 0;
	}

	int setToNeighboringBetaPair(int direction, int tid=0) {
		int steps_and_dir = 0;
		// If direction is set, it dominates
		if(direction != 0)
			steps_and_dir = direction;
		// The threads are distributed around the actual beta-values in zigzag-manor
		else {
			int remainder = tid % 2;
			steps_and_dir = tid / 2;

			if(remainder > 0)
				steps_and_dir = -steps_and_dir;
		}

		std::map<BetaPair, BKZInfo<FT1>>& the_map = this->_bench_mark->getDataMap();
		//_bench_mark->printMap();
		auto it = the_map.find(this->_bkz_info.betas);

		if(it == the_map.end()) {
			cerr << "Beta pair not registered in map. Failure." << endl;
			exit(-4200);
		}

		// Try to find an upper neighbor
		if(steps_and_dir > 0) {
			int target = steps_and_dir;
			while (target > 0 && it!=the_map.end()) {
				it++;
				target--;
			}

			// If everything was fine
			if(target==0 && it!=the_map.end()) {
				this->_bkz_info = it->second;
				//cout << "New beta pair (+): " << it->first.first << " / " << it->first.second << endl;
			}

			// Choose a random entry
			else {
				std::random_device seeder_mod;
				std::mt19937 engine(seeder_mod());
				std::uniform_int_distribution<int> dist_offset(0, the_map.size());
				int offset = dist_offset(engine);
				it=the_map.begin();
				for( int i=0; i < offset; i++) {
					it++;
				}
			}
		}

		// Try to find lower neighbor
		else {
			int target = abs(steps_and_dir);
			while (target > 0 && it!=the_map.begin()) {
				--it;
				target--;
			}

			// If everything went fine
			if(target==0) {
				this->_bkz_info = it->second;
				//cout << "New beta pair (-): " << it->first.first << " / " << it->first.second << endl;
			}

			// Choose a random entry
			else {
				std::random_device seeder_mod;
				std::mt19937 engine(seeder_mod());
				std::uniform_int_distribution<int> dist_offset(0, the_map.size());
				int offset = dist_offset(engine);
				it=the_map.begin();
				for( int i=0; i < offset; i++) {
					it++;
				}
			}
		}

		// If beta changed, the costs also changed
		this->calculateCosts();
		return 0;
	}

	FT1 setOpimtalBetaPair() {
		std::map<BetaPair, BKZInfo<FT1>>& the_map = this->_bench_mark->getDataMap();

		BetaPair bestpair;
		FT1 mincost = FT1(std::numeric_limits<long double>::max());
		vector<FT1> bestprobs;

		vector<FT1> t_reduction_tmp;
		vector<FT1>  est_nodes_tmp;
		vector<FT1>  t_enums_tmp;
		vector<FT1>  succ_prob_tmp;
		vector<FT1>  base_costs_tmp;

		t_reduction_tmp.reserve(_ainfo._number_of_random_bases);
		est_nodes_tmp.reserve(_ainfo._number_of_random_bases);
		t_enums_tmp.reserve(_ainfo._number_of_random_bases);
		succ_prob_tmp.reserve(_ainfo._number_of_random_bases);
		base_costs_tmp.reserve(_ainfo._number_of_random_bases);

		t_reduction_tmp.resize(_ainfo._number_of_random_bases);
		est_nodes_tmp.resize(_ainfo._number_of_random_bases);
		t_enums_tmp.resize(_ainfo._number_of_random_bases);
		succ_prob_tmp.resize(_ainfo._number_of_random_bases);
		base_costs_tmp.resize(_ainfo._number_of_random_bases);

		// Within each beta-pair, we have to consider all matrices contained
		for(auto it = the_map.begin(); it != the_map.end(); it++) {
			BKZInfo<FT1>& bkz_info = it->second;
			const vector<vector<FT1>>& bs_tmp = this->_bench_mark->getBsForBetaPair(bkz_info.betas);

			FT1 costs = FT1(0.0);
			for(int i=0; i<_ainfo._number_of_random_bases; i++) {
				// We can calculate bases in parallel but with a time penalty of factor 1.5
				t_reduction_tmp[i] = (FT1(bkz_info.rtime[i]));

				const long double* bs_ptr = bs_tmp[i].data();
				FT1 amax_tmp = bkz_info.amax[i];
				est_nodes_tmp[i] = estimateNumberOfNodesInTree<FT1, FT1>(_prun_func, amax_tmp, _dim, bs_ptr, _probs[i], !_probs_initialized[i]);

				_probs_initialized[i] = true;
				t_enums_tmp[i] = _ainfo._time_per_node * est_nodes_tmp[i];

				succ_prob_tmp[i] = _probs[i][_dim];

				// We do _number_of_random_bases of shots in parallel an only at least 1 has to hit
				FT1 prob_one_shot = FT1(1) - pow(FT1(1)-succ_prob_tmp[i], FT1(_ainfo._number_of_parallel_reducing_threads));

				// Like in the paper of Gama et al.
				base_costs_tmp[i] = (t_reduction_tmp[i] + t_enums_tmp[i]) / prob_one_shot;

				// With Bernoulli and 0.99 chance as average
				//_base_costs[i] = (log(FT1(0.001)) / (log(FT1(1.0) - prob_one_shot))) * (_t_reduction[i] + _t_enums[i]);

				costs += base_costs_tmp[i];
			}
			costs /= FT1(_ainfo._number_of_random_bases);

			// Store new best Solution
			if(costs < mincost && abs(costs - mincost) > 0.03 ) {
				mincost = costs;
				bestpair = it->first;
				this->_bkz_info = it->second;

				for(int i=0; i<_ainfo._number_of_random_bases; i++) {
					_succ_prob[i] = succ_prob_tmp[i];
					_est_nodes[i] = est_nodes_tmp[i];
					_t_reduction[i] = t_reduction_tmp[i];
					_t_enums[i] = t_enums_tmp[i];
					_base_costs[i] = base_costs_tmp[i];
				}
			}
		}

		this->_costs_calculated = true;
		this->_costs = mincost;

		return _costs;
	}

	FT1 calculateCosts() {
		_costs = FT1(0.0);
		for(int i=0; i<_ainfo._number_of_random_bases; i++) {
			// We can calculate bases in parallel but with a time penalty of factor 1.5
			_t_reduction[i] = (FT1(_bkz_info.rtime[i]));

			const vector<vector<FT1>>& bs_tmp = this->_bench_mark->getBsForBetaPair(_bkz_info.betas);
			const long double* bs_ptr = bs_tmp[i].data();
			FT1 amax_tmp = _bkz_info.amax[i];
			_est_nodes[i]= estimateNumberOfNodesInTree<FT1, FT1>(_prun_func, amax_tmp, _dim, bs_ptr, _probs[i], !_probs_initialized[i]);

			_probs_initialized[i] = true;
			_t_enums[i] = _ainfo._time_per_node * _est_nodes[i];

			_succ_prob[i] = _probs[i][_dim];

			// We do _number_of_random_bases of shots in parallel an only at least 1 has to hit
			FT1 prob_one_shot = _succ_prob[i] ;//FT1(1) - pow(FT1(1)-_succ_prob[i], FT1(_ainfo._number_of_instances));
			FT1 prob_thread_shots =  FT1(1) - pow((FT1(1) - prob_one_shot), _ainfo._number_of_parallel_reducing_threads);

			//cout << prob_one_shot << " vs. " << prob_thread_shots << endl;

			// Like in the paper of Gama et al.
			_base_costs[i]= (_t_reduction[i] + _t_enums[i]) / prob_thread_shots;

			// With Bernoulli and 0.99 chance as average
			//_base_costs[i] = (log(FT1(0.0001)) / (log(FT1(1.0) - prob_one_shot))) * (_t_reduction[i] + _t_enums[i]);

			_costs += _base_costs[i];
		}
		_costs /= FT1(_ainfo._number_of_random_bases);
		//cout << "Bases: " << _ainfo._number_of_random_bases << endl;
		_costs_calculated = true;
		_costs = round( _costs * 1000.0) / 1000.0;
		return this->_costs;
	}

	std::string getBetaConfigString() {
		std::string configstr = "PreBeta: ";
		configstr.append(to_string(_bkz_info.prebeta));
		configstr.append(" / Beta: ");
		configstr.append(to_string(_bkz_info.beta));

		return configstr;
	}

	int modifyToNeighborSolution(bool update_beta_pair=true) {
		std::random_device seeder_mod;
		std::mt19937 engine(seeder_mod());

		// Select the dimension to change
		// Always two neighboring values are changed at once
		std::uniform_int_distribution<int> dist_dim(0, (_dim-3) / 2);

		int change_dim = dist_dim(engine);
		// Select percentage if change Range [-1.0%, ..., +1.0%]
		std::uniform_real_distribution<double> dist_val(-0.01, 0.01);
		double perci = dist_val(engine);

		while(abs(perci) < 10e-5 || (perci < 0 && change_dim==0)) {
			perci = dist_val(engine);
		}

		//FT1 perc = (FT1(perci) / 2.0) / 100.0;
		FT1 perc = FT1(perci);

		// Make sure that the function is still monotonically increasing
		// Calculate percental difference to corresponding neighbor solution
		change_dim *=2;
		int neigh_dim = -1;
		FT1 dimval = _prun_func[ change_dim ];
		FT1 valchange;
		FT1 newval;

		if(perc > 0.0) {
			neigh_dim = change_dim + 2;
			FT1 neighval = _prun_func[ neigh_dim ];

			valchange = dimval * perc;
			newval = dimval + valchange;

			if(newval > neighval) {
				if(change_dim > 0) {
					//newval = std::max(neighval, _prun_func[change_dim-2]);
                    newval = neighval;
				}
				else
					newval=neighval;//*0.999999;
			}

		}
		else {
			neigh_dim = change_dim - 2;
			FT1 neighval = _prun_func[ neigh_dim ];

			valchange = dimval * perc;
			newval = dimval + valchange;

			if(newval < neighval) {
				//newval = std::min(neighval, _prun_func[change_dim + 2]);//*1.0000001;
                newval = neighval;
			}
		}

		_prun_func[ change_dim ] = newval;
		_prun_func[ change_dim+1 ] = newval;

		// Check feasiability
		for(long unsigned int i=0; i<_prun_func.size()-1; i++) {
				if(_prun_func[i] > _prun_func[i+1]) {
					cerr << "Prunfunc invalid! :" << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
					<< i << " (" << _prun_func[i] << " / " << _prun_func[i+1] << ")"
							<< "Perc: " << perc << " @"<< change_dim << " with neighbor "
							<< neigh_dim << " with value: " << _prun_func[ change_dim + 2 ] << endl;
					exit(-1000);
				}
		}
		//cout << "Test passed" << endl;

		// Through the change in pruning function, the probability changes
		this->_costs_calculated = false;
		for(int i=0; i<_ainfo._number_of_random_bases; i++) {
			this->_probs_initialized[i] = false;
		}

		if(update_beta_pair) {
			setOpimtalBetaPair();
		}
		else {
			calculateCosts();
		}


		return 0;
	}

protected:
	void setRandomID() {
		// Set a random id
		std::random_device seeder;
		std::mt19937 engine(seeder());
		std::uniform_int_distribution<int> dist(1, 999999);
		this->_id = dist(engine);
	}

	int _id;
	bool _costs_calculated;
	FT1 _costs;

	vector<FT1> _succ_prob;
	vector<FT1> _t_nodes;
	vector<FT1> _est_nodes;
	vector<FT1> _t_reduction;
	vector<FT1> _t_enums;
	vector<FT1> _base_costs;
	vector<FT1> _probs_initialized;
	vector<vector<FT1>> _probs;


private:
	BKZInfo<FT1> _bkz_info;
	BKZBenchmarker<FT1>* _bench_mark;
	vector<FT1> _prun_func;
	int _dim;

	AnnealInfo<FT1> _ainfo;
};



#endif /* SRC_ANNEALINGSOLUTION_HPP_ */