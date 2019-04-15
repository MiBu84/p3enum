/*
 * BKZBenchmarker.hpp
 *
 *  Created on: 25.02.2019
 *      Author: mburger
 */

#ifndef SRC_BKZBENCHMARKER_HPP_
#define SRC_BKZBENCHMARKER_HPP_

#include <map>
#include <vector>

#include <fplll.h>
#include <fplll/gso.h>
#include <NTL/mat_ZZ.h>
#include <boost/progress.hpp>
#include <omp.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace NTL;

//#define PARFACTOR 10

template <class FT>
class AnnealInfo {
public:
	int _number_of_random_bases; // Determines how many threads work in parallel on different solutions
	int _number_of_annealing_threads; // How many threads will anneal to increase number of results
	int _number_of_instances; // How many parallel BKZ-Instances will be executed in parallel afterward?
	int _number_of_enum_threads; // How many threads will afterward work on parallel bases

	FT _time_per_node; // How long does it take to calculate one node in enumeration tree

	AnnealInfo() {
		_number_of_random_bases=1;
		_number_of_annealing_threads=1;
		_number_of_instances=1;
		_number_of_enum_threads=1;
		_time_per_node=FT(-1);
	}

	AnnealInfo(const AnnealInfo& other) {
		_number_of_random_bases = other._number_of_random_bases;
		_number_of_annealing_threads = other._number_of_annealing_threads;
		_number_of_instances = other._number_of_instances;
		_number_of_enum_threads =other._number_of_enum_threads;
		_time_per_node = other._time_per_node;
	}

	AnnealInfo(AnnealInfo& other) {
		_number_of_random_bases = other._number_of_random_bases;
		_number_of_annealing_threads = other._number_of_annealing_threads;
		_number_of_instances = other._number_of_instances;
		_number_of_enum_threads =other._number_of_enum_threads;
		_time_per_node = other._time_per_node;
	}

	AnnealInfo& operator= (const AnnealInfo & other) {
		_number_of_random_bases = other._number_of_random_bases;
		_number_of_annealing_threads = other._number_of_annealing_threads;
		_number_of_instances = other._number_of_instances;
		_number_of_enum_threads =other._number_of_enum_threads;
		_time_per_node = other._time_per_node;
		return *this;
	}

	// ToDo: Is this required or covered by AnnealInfo& operator= (const AnnealInfo & other)?
	AnnealInfo& operator= (AnnealInfo & other) {
		_number_of_random_bases = other._number_of_random_bases;
		_number_of_annealing_threads = other._number_of_annealing_threads;
		_number_of_instances = other._number_of_instances;
		_number_of_enum_threads =other._number_of_enum_threads;
		_time_per_node = other._time_per_node;
		return *this;
	}
};

// <prebeta, beta>
typedef std::pair<short, short> BetaPair;

struct BetaConfig{
	short prebeta_start;
	short prebeta_end;
	short beta_start;
	short beta_end;


	// Implements the "\in" operation
	bool operator==(const BetaPair& a) const
	{
		if(
				a.first >= prebeta_start && a.first <= prebeta_end &&
				a.second >= beta_start && a.second <= beta_end)
			return true;
		return false;

	}
};





template <class FT>
struct BKZInfo{
	BetaPair betas;
	vector<double> rtime;
	vector<int> seed;
	vector<double> amax;
};

template <class FT1>
class BKZBenchmarker {

public:

	BKZBenchmarker(AnnealInfo<FT1>& ainfo) {
		this->_dim = -1;
		this->_data.clear();
		this->_ainfo = ainfo;

		this->B_org.reserve(ainfo._number_of_random_bases);
		this->B_org.resize(ainfo._number_of_random_bases);
	}
	~BKZBenchmarker() {

	}

	BKZInfo<FT1>& next() {
		auto it = _data.begin();
		return it->second;
	}

	int initialize(const mat_ZZ& B_in, const BetaConfig& bconfig) {
		ZZ_mat<mpz_t> U;

		vec_ZZ randint;
		randint.SetLength(5);
		randint[0]=-2; randint[1]=-1; randint[2]=0; randint[3]=1; randint[4]=2;
		srand (time(NULL));

		vector<mat_ZZ> Bs;
		Bs.reserve(_ainfo._number_of_random_bases);
		_dim = B_in.NumCols();
		Bs.push_back(B_in);

		// Create a bunch of randomized matrices to cover a broader range in analysis
		// Different bases must be available in inppap folder!
		int diff_bases = Configurator::getInstance().ann_num_different_bases;

		int act_seed = 0;
		double act_amax = 0.0;
		unsigned int seed_cnt = 0;

		int seed_cnt_max = (int)ceil(((double)(_ainfo._number_of_random_bases)/(double)(diff_bases)));
		mat_ZZ  B_read = B_in;
		std::string filebase = "inppap/svpc-e1-d";
		filebase.append(to_string(_dim));
		filebase.append("-");

		act_seed = Configurator::getInstance().ann_seeds_different_bases[seed_cnt];
		act_amax = Configurator::getInstance().ann_amax_different_bases[seed_cnt];

		//cout << "Number of bases: " << _ainfo._number_of_random_bases << endl;
		//cout << "Number of different bases: " << diff_bases << endl;
		//cout << "Number of bases for one seed: " << seed_cnt_max << endl;

		this->_seeds.push_back(act_seed);
		this->_amax.push_back(act_amax);
		//cout << act_seed << " / " << act_amax << endl;

		for(int i=1; i < _ainfo._number_of_random_bases;i++) {
			//cout << "i: " << i << endl;
			if(i % seed_cnt_max == 0) {

				// ToDo: Reed the correct seed from array
				seed_cnt++;

				if(seed_cnt > Configurator::getInstance().ann_seeds_different_bases.size() ||
						seed_cnt > Configurator::getInstance().ann_amax_different_bases.size()) {
					cerr << "More seeds for annealing requested than given! Sticking to seed 0" << endl;
				}

				else {
					act_seed = Configurator::getInstance().ann_seeds_different_bases[seed_cnt];
					act_amax = Configurator::getInstance().ann_amax_different_bases[seed_cnt];
				}

				//cout << "New seed: " << act_seed << endl;

				std::string filen = filebase;
				filen.append(to_string(act_seed));
				filen.append(".txt");
				std::cout << "Trying to read: " << filen << std::endl;
				if(readRandomLattice(B_read, filen) < 0)
				{
					cerr << "Reading lattice for Benchmarking failed. Exiting." << endl;
					exit(-4100);
				}
				Bs.push_back(B_read);
			}
			else {
				Bs.push_back(randomizeMatrix(B_read, randint, _dim, i));
			}
			this->_seeds.push_back(act_seed);
			this->_amax.push_back(act_amax);
			//cout << act_seed << " / " << act_amax << endl;
		}


		// Copy matrices to local fplll format
		for(int i=0; i < _ainfo._number_of_random_bases; i++) {
			B_org[i].set_cols(B_in.NumCols());
			B_org[i].set_rows(B_in.NumRows());

			for(int r=0; r < B_in.NumRows(); r++) {
				for(int c=0; c < B_in.NumCols(); c++) {
					mpz_t aa;
					mpz_init(aa);

					std::stringstream ssa;
					ssa << Bs[i][r][c];
					mpz_set_str( aa, ssa.str().c_str(),10);

					B_org[i][r][c] = aa;
				}
			}
		}


		readTablesFromFile<FT1>(bconfig);
		return doBenchmark(bconfig);
	}

	int initialize(const ZZ_mat<mpz_t> B_in, const BetaConfig& bconfig) {
		cerr << "Function initialize(const ZZ_mat<mpz_t> B_in, const BetaConfig& bconfig) not implemented" << endl;
		B_org[0] = B_in;
		this->_data.clear();
		return doBenchmark(bconfig);
	}

	int printMap() {
		for(auto it = _data.begin(); it != _data.end(); it++) {
			cout << it->second.betas.first << " , ";
			cout << it->second.betas.second << endl;

			for(int i=0; i < _ainfo._number_of_random_bases; i++) {
				//cout << it->second.prebeta << " , ";
				//cout << it->second.beta << " , ";
				cout << it->second.rtime[i] << " , ";
				cout << it->second.seed[i] << " , ";
				cout << it->second.amax[i] << endl;


				//for(auto it2 = it->second.bs[i].begin(); it2 != it->second.bs[i].end(); it2++) {
				//	cout << *it2 << endl;
				//}
			}
		}
		return 0;
	}

	const vector<vector<FT1>> getBsForBetaPair (const BetaPair& bpair) {
		return this->_bs_data[bpair];
	}

	std::string createFilenameForConfig() {
		std::string filename = "BKZPeformance_dim";
		filename.append(std::to_string(this->_dim));
		filename.append("_bases");
		filename.append(std::to_string(_ainfo._number_of_random_bases));
		if(_ainfo._number_of_random_bases != _ainfo._number_of_instances)
		{
			filename.append("_instances");
			filename.append(std::to_string(_ainfo._number_of_instances));
		}
		filename.append(".bf");
		return filename;
	}

	template <class FT>
	int writeTablesToFile() {
		ofstream myfile;
		myfile.open (createFilenameForConfig());

		for(auto it = _data.begin(); it != _data.end(); it++) {
			myfile << it->second.betas.first << " , ";
			myfile << it->second.betas.second << endl;

			for(int i=0; i < _ainfo._number_of_random_bases; i++) {
				myfile << it->second.betas.first << " , ";
				myfile << it->second.betas.second << " , ";
				myfile << it->second.rtime[i] << " , ";
				myfile << it->second.seed[i] << " , ";
				myfile << it->second.amax[i] << endl;

				const vector<vector<FT1>>& bs_temp = _bs_data[it->second.betas];
				for(auto it2 = bs_temp[i].begin(); it2 != bs_temp[i].end(); it2++) {
					myfile << *it2 << endl;
				}
			}
		}

		myfile.close();
		return 0;
	}


	template <class FT=double>
	int readTablesFromFile(const BetaConfig& bconfig) {
		string line;
		bool is_first=true;
		int number_of_base = 0;

		ifstream myfile (createFilenameForConfig());
		char delim = ',';
		std::string item;

		if (myfile.is_open()) {
			BKZInfo<FT1> tmp = BKZInfo<FT1> ();
			vector<vector<FT1>> bs_temp;

		    while ( getline (myfile,line) )
		    {
		    	std::vector<std::string> elems;
		    	std::stringstream ss(line);

				while (std::getline(ss, item, delim)) {
					elems.push_back(std::move(item)); // if C++11
				}

				if(elems.size() == 1) {
					bs_temp[number_of_base].push_back(FT(std::stod(elems[0])));
				}

				// New random base with same beta-config
				else if(elems.size() >= 3) {
					number_of_base++;
					tmp.rtime[number_of_base] = std::stod(elems[2]);

					// Read the Amax for the corresponding dimension AND seed
					if(elems.size() == 5) {
						tmp.seed[number_of_base] = std::stoi(elems[3]);
						tmp.amax[number_of_base] = std::stod(elems[4]);
					}
					// For compatibiliy before 2019-03-22
					else {
						tmp.seed[number_of_base] = 0;
						tmp.amax[number_of_base] = Configurator::getInstance().Amax;
					}
				}

				// Save old bkz Info and create new one
				else if(elems.size() == 2) {
					if(!is_first)
					{
						BetaPair idx = BetaPair(tmp.betas.first , tmp.betas.second);

						if(bconfig == idx)
						{
							this->_data[idx] = tmp;
							this->_bs_data[idx] = bs_temp;
						}
					}

					// New entry;
					tmp = BKZInfo<FT1> ();

					tmp.rtime.reserve(_ainfo._number_of_random_bases);
					tmp.rtime.resize(_ainfo._number_of_random_bases);
					tmp.seed.reserve(_ainfo._number_of_random_bases);
					tmp.seed.resize(_ainfo._number_of_random_bases);
					tmp.amax.reserve(_ainfo._number_of_random_bases);
					tmp.amax.resize(_ainfo._number_of_random_bases);
					bs_temp.clear();
					bs_temp.reserve(_ainfo._number_of_random_bases);
					bs_temp.resize(_ainfo._number_of_random_bases);

					tmp.betas.first = std::stoi(elems[0]);
					tmp.betas.second = std::stoi(elems[1]);
					number_of_base = -1;
					is_first = false;
				}

				else {
					cerr << "Wrong length in data map file" << endl;
				}
		    }

		    // Insert the final entry
			BetaPair idx = BetaPair(tmp.betas.first , tmp.betas.second);

			// Only store if beta-pair lies in the range of the BKZInfo
			if(bconfig == idx)
			{
				this->_data[idx] = tmp;
				this->_bs_data[idx] = bs_temp;
			}

			cout << "Successfully read " << _data.size() << " entries from file." << endl;
		    myfile.close();
		    return 0;
		}

		else {
			cout << "Unable to open map file";
			return -3333;
		}
	}

	int doBenchmark(const BetaConfig& bconfig) {
		double benchstart = omp_get_wtime();
		cout << "Benchmarking with " << _ainfo._number_of_instances << " threads." << endl;
		boost::progress_display show_progress(
				(((bconfig.prebeta_end - bconfig.prebeta_start)/2 + 1) * ((bconfig.beta_end - bconfig.beta_start)/2 + 1))*_ainfo._number_of_random_bases);

		int already_read = 0;
		int calculated_new = 0;
		vector<Strategy> strategies_full;
		strategies_full = load_strategies_json(strategy_full_path("default.json"));

		for(short prebet = bconfig.prebeta_start; prebet <= bconfig.prebeta_end; prebet+=2) {
			for(short bet = bconfig.beta_start; bet <=bconfig.beta_end; bet+=2) {
				BKZInfo<FT1> binfo;
				binfo.betas.first = prebet;
				binfo.betas.second = bet;
				binfo.rtime.reserve(_ainfo._number_of_random_bases);
				binfo.rtime.resize(_ainfo._number_of_random_bases);
				binfo.seed.reserve(_ainfo._number_of_random_bases);
				binfo.seed.resize(_ainfo._number_of_random_bases);
				binfo.amax.reserve(_ainfo._number_of_random_bases);
				binfo.amax.resize(_ainfo._number_of_random_bases);
				vector<vector<FT1>> bs_temp;
				bs_temp.reserve(_ainfo._number_of_random_bases);
				bs_temp.resize(_ainfo._number_of_random_bases);
				std::vector<ZZ_mat<mpz_t> > reduced_bases;
				reduced_bases.clear();
				reduced_bases.reserve(_ainfo._number_of_random_bases);
				reduced_bases.resize(_ainfo._number_of_random_bases);
				vector<double> tim;
				tim.reserve(_ainfo._number_of_random_bases);
				tim.resize(_ainfo._number_of_random_bases);

				// Check if already read, then it needs not to be produced again
				BetaPair testpair = BetaPair(prebet, bet);
				if(this->_data.count ( testpair ) ) {
					already_read++;

					for(int i=0; i < _ainfo._number_of_random_bases; i++)
						++show_progress;
					continue;
				}

		#pragma omp parallel for num_threads( _ainfo._number_of_instances )
				for(int i=0; i < _ainfo._number_of_random_bases; i++) {
					// Make temporary copy of input matrix
					ZZ_mat<mpz_t> B_loc = B_org[i];
					ZZ_mat<mpz_t> empty1, empty2;

					double tstart = omp_get_wtime();
					vector<Strategy> strategies;
					if(prebet > 0)
					{
						// BKZ 1.0 reducion
						BKZParam param20(prebet, strategies);
						if(prebet > 20) {
							param20.flags |=  BKZ_AUTO_ABORT;
						}
						param20.delta = Configurator::getInstance().glob_delta;
						bkz_reduction 	(&B_loc, &empty1,  param20, FT_DEFAULT, 0);
						//LLLReduction<Z_NR<mpz_t>, FP_NR<FT1>> lll_obj(mygso2, param20.delta, LLL_DEF_ETA, LLL_DEFAULT);
						//BKZReduction<Z_NR<mpz_t>, FP_NR<FT1>> bkz_obj1(mygso2, lll_obj, param20);
						//bkz_obj1.bkz();
					}

					// BKZ 2.0 reduction
					BKZParam paramval(bet, strategies_full);
					paramval.delta = Configurator::getInstance().glob_delta;
					paramval.flags |= BKZ_AUTO_ABORT;
					bkz_reduction 	(&B_loc, &empty2,  paramval, FT_DEFAULT, 0);
					//BKZReduction<Z_NR<mpz_t>, FP_NR<FT1>> bkz_obj2(mygso2, lll_obj, paramval);
					//bkz_obj2.bkz();

					double tend = omp_get_wtime();

#pragma omp critical
{
					//double regstart = omp_get_wtime();
					++show_progress;
					reduced_bases[i] = B_loc;

					tim[i] = tend - tstart;
}
				} // PARFACTOR

				// Write stuff in serial, perhaps this works without crash
				for(int i=0; i < _ainfo._number_of_random_bases; i++) {
					ZZ_mat<mpz_t> empty3, empty4;
					MatGSO<Z_NR<mpz_t>, FP_NR<long double>> mygso(reduced_bases[i], empty3, empty4, 0);
					mygso.update_gso();

					// Copy bvals for storage
					FP_NR<long double> tmp;

					bs_temp[i].reserve(reduced_bases[i].get_rows());
					binfo.rtime[i] = tim[i];
					binfo.seed[i] = this->_seeds[i];
					binfo.amax[i] = this->_amax[i];
					//cout << "Seed " << binfo.seed << " / Amax " << binfo.amax << endl;
					for(int ii=0; ii < reduced_bases[i].get_rows(); ii++) {
						tmp = mygso.get_r(tmp, ii, ii);
						bs_temp[i].push_back(FT1(tmp.get_data()));
					}
				}

				BetaPair npair(prebet, bet);
				_data[npair] = binfo;
				_bs_data[npair] = bs_temp;
				calculated_new++;
			}
		}
		cout << endl << "Size of map: " << _data.size() << " (Calc: " << calculated_new
				<< " / Read: " << already_read << ")"
				<< " in " <<  omp_get_wtime() - benchstart << " seconds." << endl;
		//printMap();
		writeTablesToFile<double>();
		return 0;
	}

	std::map<BetaPair, BKZInfo<FT1>>& getDataMap() {
		return this->_data;
	}

	unsigned int size() {
		return this->_data.size();
	}

	int getSeedFor(int pos) {
		return _seeds[pos];
	}

	double getAmaxFor (int pos) {
		return _amax[pos];
	}

private:
	std::map<BetaPair, BKZInfo<FT1>> _data;
	std::map<BetaPair, vector<vector<FT1>>> _bs_data;

	vector<ZZ_mat<mpz_t>> B_org;
	int _dim;
	AnnealInfo<FT1> _ainfo;
	vector<int> _seeds;
	vector<double> _amax;
};



#endif /* SRC_BKZBENCHMARKER_HPP_ */
