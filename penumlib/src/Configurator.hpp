/*
 * Configurator.hpp
 *
 *  Created on: 12.07.2018
 *      Author: mburger
 */

#ifndef SRC_CONFIGURATOR_HPP_
#define SRC_CONFIGURATOR_HPP_

#include <cstring>
#include <string>
#include <limits>
#include <omp.h>

#include "pruningfunc.hpp"

class Configurator {
public:
	static Configurator& getInstance()
	{
		static Configurator instance;
		return instance;
	}

	static Configurator* _instance;

	// Configuration variables
	double maxas[200];

	// Basis reduction
	bool dolll;
	string reducedfilepath;

	int forced_height;
	double glob_delta;
    int BKZinstances; // How many threads are employed for parallel instances of BKZ
    int Enumthreads; // Limits number of threads in enumeration. Default: All.
	int glob_beta;
    int prebeta; // How high is the beta for a pre-run of BKZ without pruning, but with auto abort
    
	bool use_ntl_enum;
	int par_threshold;
	int cand_queue_size;
	int bkz_height;
	double gauss_estimate;


	bool auto_height_tuning;
	bool do_gauss_estimate;
	bool do_bkz_auto_tuning;
	bool use_fplll_bkz;

	Pruning enum_prune;
	double prune_param;
	vector<double> ext_pruning_function;
	string ext_pruning_function_file;

	double Amax;
	int trials; // How many trial of pruned and randomized repetitions

	double Aend;
	bool iterative_enumeration;

private:
	Configurator() {
		dolll=true;
		reducedfilepath="";

		forced_height = -1;
		glob_delta = 0.9;
		glob_beta = -1;
        prebeta = 18;        
        
		use_ntl_enum = true;
		par_threshold = 10;
		cand_queue_size = 250;
		auto_height_tuning = true;
		do_gauss_estimate = true;
		gauss_estimate = 0.0;
		do_bkz_auto_tuning = false;
		use_fplll_bkz = true;
		bkz_height = 10;

		enum_prune = NO;
		prune_param = 1.00;
		Amax = numeric_limits<double>::max();
		trials = 40; // Following Paper of Schneider and Dagdalen this garantuees success of 99%
        
        #pragma omp parallel
        #pragma omp single
        {
            Enumthreads = omp_get_num_threads();
            BKZinstances = -1; // Causes program to execute half as many BKZ instances as threads configured for enum
        }

		Aend = numeric_limits<double>::max();
		iterative_enumeration = false;

		// Default: Use scaled pruning function
		ext_pruning_function_file = "NOT";

		// Setting max A's from the svpchallenge website
		maxas[50] = 1881;
		maxas[51] = 1885;
		maxas[52] = 1906;
		maxas[53] = 1896;
		maxas[54] = 1921;
		maxas[55] = 1973;
		maxas[56] = 1974;
		maxas[57] = 1953;
		maxas[58] = 2017;
		maxas[59] = 2069;
		maxas[60] = 1943;
		maxas[61] = 2112;
		maxas[62] = 2092;
		maxas[63] = 2056;
		maxas[64] = 2103;
		maxas[65] = 2095;
		maxas[66] = 2099;
		maxas[67] = 2187;
		maxas[68] = 2141;
		maxas[69] = 2212;
		maxas[70] = 2143;
		maxas[71] = 2228;
		maxas[72] = 2179;
		maxas[73] = 2190;
		maxas[74] = 2254;
		maxas[75] = 2280;
		maxas[76] = 2241;
		maxas[77] = 2229;
		maxas[78] = 2273;
		maxas[79] = 2346;
		maxas[80] = 2272;
		maxas[82] = 2303;
		maxas[83] = 2304; // seed=64

		maxas[85] = 2347; // seed=64
		maxas[86] = 2357;

		maxas[93] = 2396; // seed=64
		maxas[95] = 2523; // Artur
	}

	Configurator(const Configurator&);
	~Configurator() {}



};

#endif /* SRC_CONFIGURATOR_HPP_ */
