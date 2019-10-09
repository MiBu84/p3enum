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
#include <boost/algorithm/string/classification.hpp> // Include boost::for is_any_of
#include <boost/algorithm/string/split.hpp>

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

	int readSolutionFromFile(string filepath, const int dim, const int bkz_entry) {
		//cout << "Called Reader from " << bkz_entry << endl;
		std::string locpath = filepath;
		std::ifstream file(locpath.c_str());

		if (!file)
		{
		    cerr << "No valid pruningfile to continue from found " << locpath << endl;
		    return 11001;
		}

		std::string   line;

		string str1="PreBeta:";
		string str2="final";

		int entrycount = 0;
		int act_prun_entry = -1;
		int mode = -1;

		while(std::getline(file, line))
		{
			//cout << line << endl;
			std::vector<std::string> words;
			boost::split(words, line, boost::is_any_of(", "), boost::token_compress_on);

			if(words.size() == 0)	continue;

			// The last solution was processed, now only the best solution follows
			if(words.size() > 1) {
				if(str2.compare(words[1]) == 0) {
						break;
					}
			}

			if (str1.compare(words[0]) == 0) {
				mode = 1;
			}

			// Check if this is the entry, we are interested in
			// If not, skip to next entry
			if(mode==1 && entrycount != bkz_entry) {
				mode = -1;
				entrycount++;
			}

			// Check betas for correctness
			if(mode==1) {
				cout << "I found a solution for instance number " << entrycount << endl;
				if(words.size() != 5) {
					cerr << "Configfile Beta-line invalid" << endl;
					return 11002;
				}

				int prebeta = (int) std::stoi(words[1]);
				int beta = (int) std::stoi(words[4]);

				if(prebeta != this->_bkz_info.betas.first || beta != this->_bkz_info.betas.second) {
					cerr << "Betas differ: " << prebeta << " / " << this->_bkz_info.betas.first
							<< " ; " << beta << " / " << this->_bkz_info.betas.second << endl;
					return 11003;
				}
				mode = 2;
			}
			// The next line gives the BKZ Times, should be read from the BKZ Benchmarker
			else if (mode == 2) {
				//cout << words[0] << " " << words[1] << endl;
				mode = 3;
			}

			// The next line gives the number of nodes, should be recalculated below
			else if (mode == 3) {
				//cout << words[0] << " " << words[1] << endl;
				mode = 4;
			}

			// The next line gives the time of ENUM, should be recalculated below
			else if (mode == 4) {
				//cout << words[0] << " " << words[1] << endl;
				mode = 5;
			}

			// The next line gives the success probs, should be recalculated below
			else if (mode == 5) {
				//cout << words[0] << " " << words[1] << endl;
				mode = 6;
			}

			// The next line gives the total costs, should be recalculated below
			else if (mode == 6) {
				//cout << words[0] << " " << words[1] << endl;
				mode = 7;
				this->_prun_func.clear();
				this->_prun_func.reserve(dim);
				this->_prun_func.resize(dim);
				act_prun_entry = dim - 1;
			}

			// Now the actual pruning function starts
			else if (mode==7) {
				if(words.size() != 1) {
					cerr << "Error in prunfunc in file" << endl;
					return 11004;
				}
				this->_prun_func[act_prun_entry] = (FT1)std::atof(words[0].c_str());
				//cout << words[0].c_str() << endl;
				act_prun_entry--;

				if(act_prun_entry < 0) {
					mode=8;
				}
			}

			else if (mode==8) {
				for (auto it = _probs_initialized.begin(); it != _probs_initialized.end(); ++it) {
					*it = false;
				}

				this->calculateCosts();
				cout << "Costs are now: " << this->getCost()
						<< " before: " << words[2] <<
						endl;
				mode = -1;
				entrycount++;
			}
		}

		return 0;

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

		ss << "Overall cost: " << _costs << endl;
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
			// Uneven
            if(length%2 != 0) {
                if(ass_num == 0) {
                	ass_num = middle;
                }
                
                else if(ass_num % 2 != 0) {
                    ass_num = middle - (ass_num+1) / 2;
                }
                
                else  {
                    ass_num = middle + ass_num/2;
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

			// If each thread enumerates for itself
			bool parstrat_instance = (Configurator::getInstance().par_threshold > _dim);
			if (parstrat_instance) {
				FT1 prob_one_shot = _succ_prob[i];
				FT1 prob_thread_shots =  FT1(1) - pow((FT1(1) - prob_one_shot), _ainfo._number_of_parallel_reducing_threads);
				
				// New trial including parallelism
				// How many rounds do we expect in average?
				int expec_rounds = std::max<int>((int)((log(FT1(0.01)) / (log(FT1(1.0) - prob_thread_shots)))),1);
				_base_costs[i] = ((_t_reduction[i] + _t_enums[i]) * (double)expec_rounds);

				// Like in the paper of Gama et al.
				//_base_costs[i]= (_t_reduction[i] + _t_enums[i]) / prob_thread_shots;

				// With Bernoulli and 0.99 chance as average
				//_base_costs[i] = (log(FT1(0.0001)) / (log(FT1(1.0) - prob_thread_shots))) * (_t_reduction[i] + _t_enums[i]);

			}

			// n-parallel BKZs and then parallelized enumerations one after the other
			else {
				// How many trials do we expect?
				int expec_trials = (log(FT1(0.0001)) / (log(FT1(1.0) - _succ_prob[i])));

				// Since several BKZ can be done in parallel, less calls compared to enumeration will be required
				int expec_bkzs = expec_trials / _ainfo._number_of_parallel_reducing_threads;
				_base_costs[i] = expec_bkzs * _t_reduction[i] + expec_trials * _t_enums[i];

				//cout << "Expected Enums: " << expec_trials << " / expected BKZs: " << expec_bkzs
				//		<< " / expected time: " << _base_costs[i] << "s." << endl;
			}

			_costs += _base_costs[i];
		}
		_costs /= FT1(_ainfo._number_of_random_bases);
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
		std::uniform_real_distribution<double> dist_val(-0.1, 0.1);
		double perci = dist_val(engine);

		while(abs(perci) < 2e-2 || (perci < 0 && change_dim==0)) {
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

	BetaPair& getBetaPair() {
		return this->_bkz_info.betas;
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
