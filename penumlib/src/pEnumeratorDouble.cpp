/*
 * pEnumeratorDouble.cpp
 *
 *  Created on: 22.06.2018
 *      Author: mburger
 */


#include "pEnumeratorDouble.hpp"
#include "Configurator.hpp"
#include "Utils.hpp"
#include "MBVec.hpp"
#include "extremePruningFunction.hpp"
#include "pruningfunc.hpp"
#include "VectorStorage.hpp"
#include "SimulatedAnnealer.hpp"

#include <iostream>
#include <iomanip>      // std::setprecision
#include <math.h>
#include <omp.h>
#include <limits>
#include <fstream>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <ctime>
#include <vector>

// fpllltest
#include <cstdlib>
#include <fplll.h>

#define auxCeil(x) ((double)(long)((x)+1))

void printVector(double* v, int s, int e) {
	for(int i=s; i<=e; i++)
	cout	 << v[i] << " ";
	cout << endl;
}


pEnumeratorDouble::pEnumeratorDouble() {
	cout << "Constructor of pEnumeratorDouble called without dimension parameter!"
			<< endl;
	cand_umin=NULL;
	cand_v=NULL; // next coordinate
	cand_l=NULL; // actual costs
	cand_c=NULL; // centers of all levels in tree
	cand_Delta=NULL; // sign
	cand_delta=NULL; // offset
	cand_s = -1; // highest non-zero entry
	cand_t = -1; // actual visited level in tree
	cand_r = NULL;
	cand_sigma = NULL;


	umin = NULL;
	v = NULL;
	l = NULL;
	c = NULL;

	// To decide the zigzag pattern
	Delta  = NULL;
	delta = NULL;

	_dim = -1;
	_cnt = 2;
	_cnt2 = 0;
    _cand_cnt = 0;
	_num_threads = 1;

	candidates_queue = NULL;
	tstart_bench = NULL;
	node_cnt = 0;

	r = NULL;
	//sigma.resize(boost::extents[2][2][2]);
	sigma = NULL;
}

pEnumeratorDouble::pEnumeratorDouble(int dim) {
	_dim = dim;
	_num_threads = omp_get_max_threads();

	cand_umin=new int[dim + 3];
	cand_v=new int[dim + 3]; // next coordinate
	cand_l=new double[dim + 3]; // actual costs
	cand_c=new double[dim + 3]; // centers of all levels in tree
	cand_Delta=new int[dim + 3]; // sign
	cand_delta=new int[dim + 3]; // offset
	cand_s = 0; // highest non-zero entry
	cand_t = 0; // actual visited level in tree

	cand_r = new int[dim + 3];
	cand_sigma = new double*[dim + 3];
	for(int i=0; i<dim+3; i++) {
		cand_sigma[i] = new double[dim+3];
	}

	for(int i=0; i< _dim + 3; i++) {
		cand_r[i] = i;
		for(int j=0; j < _dim + 3; j++) {
			cand_sigma[i][j] = 0.0;
		}
	}

	//candidates_queue = MBVecQueue2<double>(Configurator::getInstance().cand_queue_size);
	candidates_queue = new MBVecQueue3<int>(Configurator::getInstance().cand_queue_size);// storage for the candidates for SVs

	prunfunc.SetLength(_dim+3);
	_conf = pruning_conf {NO, 1.05};

	for(int i=0; i<=_dim+1; i++) {
		cand_umin[i] = 0;
		cand_v[i] = 0;
		cand_l[i] = 0;
		cand_c[i] = 0;
		cand_Delta[i] = 0;
		cand_delta[i] = 1;

		prunfunc[i] = 0;
	}

	umin = new int*[_num_threads + 3];
	v = new int*[_num_threads + 3];
	l = new double*[_num_threads + 3];
	c = new double*[_num_threads + 3];
	Delta = new int*[_num_threads + 3];
	delta = new int*[_num_threads + 3];

	for(int i=0; i < _num_threads+2; i++) {
		umin[i] = new int[_dim + 3];
		v[i] = new int[_dim + 3];
		l[i] = new double[_dim + 3];
		c[i] = new double[_dim + 3];
		Delta[i] = new int[_dim + 3];
		delta[i]  = new int[_dim + 3];

		for(int j=0; j <_dim+2; j++) {
			umin[i][j] = 0;
			v[i][j] = 0;
			l[i][j] = 0.0;
			c[i][j] = 0.0;
			Delta[i][j] = 0;
			delta[i][j] = 0;
		}
	}

	// Values for center-caching within enumeration
	//this->r = new int*[this->_num_threads];
	//this->sigma = new double[this->_num_threads][_dim + 4][dim + 4];
	//this->sigma.resize(boost::extents[this->_num_threads][_dim + 4][dim + 4]);

	//for (int nt = 0; nt < this->_num_threads; nt++) {
	//	this->r[nt] = new int[ _dim + 4];
	//}

	this->r = new int*[this->_num_threads];
	this->sigma = new double**[this->_num_threads];

	for (int nt = 0; nt < this->_num_threads; nt++) {
		this->r[nt] = new int[ _dim + 4];
		this->sigma[nt] = new double*[ _dim + 4 ];
		for(int i = 0; i < _dim + 4; i++) {
			this->sigma[nt][i] = new double[ _dim + 4 ];
		}
	}

	_cnt = 1;
	_cnt2 = 0;
    _cand_cnt = 0;
	tstart_bench = NULL;
	node_cnt = NULL;
}

pEnumeratorDouble::~pEnumeratorDouble() {
	//cout << "Calling pEnumeratorDouble destructor" << endl;
	if(_dim > 0)
	{
		delete[] cand_umin;
		delete[] cand_v;
		delete[] cand_l;
		delete[] cand_c;
		delete[] cand_Delta;
		delete[] cand_delta;

		for(int i=0; i < _num_threads+2; i++) {
			delete[] umin[i];
			delete[] v[i];
			delete[] l[i];
			delete[] c[i];
			delete[] Delta[i];
			delete[] delta[i];
		}

		delete[] umin;
		delete[] v;
		delete[] l;
		delete[] c;
		delete[] Delta;
		delete[] delta;
		delete candidates_queue;



		// Values for center-caching within enumeration
		/*if(this->r != NULL) {
			for(int nt=0; nt < this->_num_threads; nt++ ) {
				delete[] this->r[nt];
				for(int i = 0; i < this->_dim + 4; i++) {
					delete[] this->sigma[nt][i];
				}
				delete[] this->sigma[nt];
			}
			delete[] this->r;
			delete[] this->sigma;
		}*/
	}
}


double pEnumeratorDouble::solveSVPMP(mat_ZZ& B, vec_ZZ& vec) {
	// Annealing mode was chosen
	if(Configurator::getInstance().do_annealing) {
		AnnealInfo<long double> ainfo;
		ainfo._number_of_random_bases = Configurator::getInstance().ann_bkz_instances;
		ainfo._number_of_annealing_threads = Configurator::getInstance().ann_annealing_threads;
		ainfo._number_of_parallel_reducing_threads = Configurator::getInstance().ann_parallel_reducing_threads;
		ainfo._number_of_different_bases = Configurator::getInstance().ann_num_different_bases;
		ainfo._time_per_node = Configurator::getInstance().ann_time_per_node; //4.145655e-8; //  3.082812277e-9

		SimulatedAnnealer<long double> annealer = SimulatedAnnealer<long double>(ainfo);
		annealer.anneal(B);
	}

		// Use external pruning functions
		if (Configurator::getInstance().ext_pruning_function_file.compare("NOT") != 0) {
			cout << "Reading pruning function from " << Configurator::getInstance().ext_pruning_function_file << endl;
			readPruning(Configurator::getInstance().ext_pruning_function_file);
		}

		vec_ZZ randint;
		randint.SetLength(5);
		randint[0]=-2; randint[1]=-1; randint[2]=0; randint[3]=1; randint[4]=2;

		double tstart_oall = omp_get_wtime();
		std::vector<mat_ZZ> Bcands;
		//std::vector<int> Bcandsorder;

		int numt = 1;
		int procnum = 0;
		bool iterative_enumeration = Configurator::getInstance().iterative_enumeration;
        double BKZmaxtime = 9000;
        bool shortest_time_set = false;

		int dimen = B.NumCols();
		double target_a = -1;
		double act_A = -1;

		double theAval = Configurator::getInstance().Amax;
		if(theAval != numeric_limits<double>::max() && theAval > 0)
			target_a = theAval * theAval;
		else if (Configurator::getInstance().do_gauss_estimate && Configurator::getInstance().gauss_estimate > 1e-5) {
			target_a = Configurator::getInstance().gauss_estimate;
		}
        else if(theAval < 0) {
           target_a = calcGaussHeuristicChallenge<ZZ>(B, dimen);
           target_a = target_a * target_a;
        }

        act_A = target_a;

    	_conf = pruning_conf {Configurator::getInstance().enum_prune,
    		Configurator::getInstance().prune_param};
    	updatePruningFuncLoc(prunfunc.data(), _conf, act_A, dimen, 0, dimen-1);

#pragma omp parallel
#pragma omp single
		{
			numt = Configurator::getInstance().BKZinstances;
            if(numt<1) {
                numt = omp_get_num_threads();
                if(numt<1) numt=1;
            }

            else {
                numt= min(omp_get_num_threads(), Configurator::getInstance().BKZinstances);
            }


			if(numt<1) numt=1;

			cout << "Preprocessing and ENUM with target length "
					<< sqrt(act_A) << " and "
					<< numt << " threads." << endl;

			Bcands.reserve(numt);
			Bcands.resize(numt);
			//Bcandsorder.reserve(numt);
			//Bcandsorder.resize(numt);
		}
		srand (time(NULL));

        omp_set_num_threads(numt); // Use pre-defined threads for all consecutive parallel regions
        int do_searching = 1;
        int act_trial = 0;

#pragma omp parallel num_threads(numt)
{
        int tid = omp_get_thread_num();
        double act_A_t = act_A;

		mat_RR mu;
		vec_RR bstar;

		double** p_mu = new double*[B.NumRows()];
		for(int i=0; i<B.NumRows(); i++)
			p_mu[i] = new double[B.NumRows()];

		double* p_bstar = new double[B.NumRows()];
		int* p_u = new int[B.NumRows()+2];
		double* p_prunfunc = new double[B.NumRows()+1];

        while(do_searching > 0)
	    //for(int trial=0; trial < Configurator::getInstance().trials; trial++)
	    {
#pragma omp atomic
        	act_trial++;

			// Reset random basis
			Bcands[tid] = B;
			//Bcandsorder[tid] = tid;
			if(omp_get_thread_num() == 0 && ((act_trial) / (numt+1))  < 1 && Configurator::getInstance().test_original_basis) {
				// Keep basis unrandomized
				cout << "Keeping org. base by thread " << omp_get_thread_num()  << endl << std::flush;
			}
			else	{
				int seed = ((act_trial+1)*numt + (tid+1)*time(NULL));
				Bcands[tid] = randomizeMatrix(Bcands[tid], randint, Bcands[tid].NumCols(),seed);
			}

			// Do BKZ in parallel and generate a bunch
			int prebeta=Configurator::getInstance().prebeta;
			if(Configurator::getInstance().dolll)
			{
				for(int bval=Configurator::getInstance().glob_beta; bval<=Configurator::getInstance().glob_beta; bval++) {
					double tstartbkz = omp_get_wtime();
					double tbkzexact_start = 1;
					double tbkzexact_end = 0;
					double tbkzheu_start = 1;
					double tbkzheu_end = 0;

					// Alternative calculation in fpllls lib
					if(Configurator::getInstance().use_fplll_bkz)
					{
						// RAM based
						ZZ_mat<mpz_t> m_tored;
						ZZ_mat<mpz_t> U;

						m_tored.set_cols(B.NumCols());
						m_tored.set_rows(B.NumRows());

						for(int r=0; r < B.NumRows(); r++) {
							for(int c=0; c < B.NumCols(); c++) {
								mpz_t aa;
								mpz_init(aa);

								std::stringstream ssa;
								ssa << Bcands[tid][r][c];
								mpz_set_str( aa, ssa.str().c_str(),10);

								m_tored[r][c] = aa;
							}
						}

						vector<Strategy> strategies;


						// First run of pure BKZ with small beta
						double bkzexacttime_start = omp_get_wtime();
						if(Configurator::getInstance().prebeta > 0)
						{
							BKZParam param20(prebeta, strategies);
							if(prebeta > 20) {
								//param20.flags |=  BKZ_MAX_TIME;// | BKZ_AUTO_ABORT; BKZ_VERBOSE;// |
								param20.flags |=  BKZ_AUTO_ABORT;
							}
							param20.delta = Configurator::getInstance().glob_delta;
							//param20.max_loops = 1;
							//param20.max_time = BKZmaxtime;
							tbkzexact_start = omp_get_wtime();
							bkz_reduction 	(&m_tored, &U,  param20, FT_DEFAULT, 0);
							tbkzexact_end = omp_get_wtime();
						}
						double bkzexacttime_end = omp_get_wtime();

						if(!shortest_time_set)
						{
							shortest_time_set = true;
							BKZmaxtime = 2.1 * (bkzexacttime_end - bkzexacttime_start);
							cout << "Setting BKZ-Max-Time to " << BKZmaxtime << " s." << endl;
						}

						// Second run with BKZ 2.0
						tbkzheu_start = omp_get_wtime();
						if(Configurator::getInstance().glob_beta > 0)
						{
							strategies = load_strategies_json(strategy_full_path("default.json"));
							BKZParam paramval(bval, strategies);
							paramval.delta = Configurator::getInstance().glob_delta;
							paramval.flags |= BKZ_AUTO_ABORT;

							//for(int beta_i=4; beta_i<=Configurator::getInstance().glob_beta; beta_i+=4)
							{
								//paramval.block_size = beta_i;
								//bkz_reduction 	(&m_tored, &U,  paramval, FT_DEFAULT, 0);
								bkz_reduction 	(&m_tored, &U,  paramval, FT_DEFAULT, 0);
							}
						}
						tbkzheu_end = omp_get_wtime();

						for(int r=0; r < B.NumRows(); r++) {
							for(int c=0; c < B.NumCols(); c++) {
								std::stringstream ssb;
								ssb << m_tored[r][c];
								Bcands[tid][r][c] = conv<ZZ>(ssb.str().c_str());
							}
						}
					}
					double tendbkz = omp_get_wtime();


					std::cout << "[Thread: " << tid << "]\t"
							<< "|b_0| of BKZ-" << prebeta << "/" << bval  << " reduced basis "
							<< lengthOfVector<ZZ, RR> (Bcands[tid][0])
							<< "(" << tendbkz-tstartbkz << "s "
							<< "[" << tbkzexact_end - tbkzexact_start << " s / "
							<<  tbkzheu_end - tbkzheu_start << "s ]"
							<< ")."
							<< std::endl;
				}
			} // Todo: Should be the same until here like in solveSVP


			// If wished for debug, the first basis is replaced by an external reduced one
			if(tid==0) {
				if(Configurator::getInstance().prereduced_base != "") {
					//readRandomLattice(Bcands[0], Configurator::getInstance().prereduced_base);
				}
			}

			// Now each thread processes the base it has reduced

			// Calculate Gram-Schmidt
			GSO(Bcands[tid], mu, bstar);

			for(int i=0; i < B.NumCols(); i++) {
				p_u[i] = 0;
				p_bstar[i] = 0;
				p_prunfunc[i] = 0;
			}

			for(int i = 0; i < B.NumRows(); i++) {
				p_bstar[i] = NTL::conv<double>(bstar[i]);
				for(int j = 0; j < B.NumCols(); j++) {
					p_mu[i][j] = NTL::conv<double>(mu[i][j]);
				}
			}
			int dim = B.NumCols();

			act_A_t = EnumDouble(p_mu, p_bstar, p_u, 0, dim-1, dim, act_A_t, 0);
			procnum++;
			NTL::Mat<long int> Bd; Bd.SetDims(B.NumRows(), B.NumCols());

			for(int i = 0; i < B.NumRows(); i++) {
				for(int j = 0; j < B.NumCols(); j++) {
					Bd[i][j] = NTL::conv<long int>(Bcands[tid][i][j]);
				}
			}

			// Standard behavior of Extreme Pruning
			// End if a vector smaller than target found
			if(act_A_t < target_a && !iterative_enumeration) {
				double tend_oall = omp_get_wtime();
				cout << "Finishing after " << act_trial+1 << " attemps in "
						<< tend_oall - tstart_oall << " s (processed "
						<< procnum << " bases)."
						<< endl;

				NTL::Vec<double> sol;	sol.SetLength(B.NumCols());
				for(int ii=0; ii<B.NumCols(); ii++)
					sol[ii] = p_u[ii];

				calcVectorLengthDouble<long int, double>(sol, Bd);
				act_A = act_A_t;
				do_searching = 0;

		    	if(Configurator::getInstance().output_success_base) {
		    		string filename = "";
		    		if(Configurator::getInstance().prereduced_base=="") {
						auto t = std::time(nullptr);
						auto tm = *std::localtime(&t);
						std::ostringstream oss;
						oss << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
						auto str = oss.str();

						filename = "Dimension-" + std::to_string(dim) + "-A-" + to_string((long)(act_A))
								+ "-" + str + ".txt";
		    		}
		    		else {
		    			filename = Configurator::getInstance().prereduced_base;
		    		}
		    		writeNTLLatticeToFile(Bcands[tid], filename);
		    	}

				exit(5);

			}

			// Iterative behavior: Search as long as there is no better solution
			// than minimum
			else if(iterative_enumeration && act_A < target_a && do_searching > 0) {
				// Also the target value given by user was found
				if(act_A_t < Configurator::getInstance().Aend * Configurator::getInstance().Aend) {
					double tend_oall = omp_get_wtime();
					cout << "Finishing iterative search after " << act_trial+1 << " attemps in "
							<< tend_oall - tstart_oall << " s (processed "
							<< procnum << " bases)."
							<< endl;
					act_A = act_A_t;
					do_searching = 0;
				}

				// The search has to go on. A short vector was found, but not the target
				act_trial = 0; // Reset the trial count
			}

			if(act_trial >= Configurator::getInstance().trials && do_searching > 0)
			{
				do_searching = -1;
			    double tend_oall = omp_get_wtime();
			    	cout << "Unsuccessful after " << Configurator::getInstance().trials << " attemps in "
			    			<< tend_oall - tstart_oall << " s (processed "
							<< procnum << " bases)."
			    			<< endl;
			    	act_A = -1;
			}
	    }
		delete[] p_bstar;
		delete[] p_u;
		delete[] p_prunfunc;
		for(int i=0; i<B.NumRows(); i++)
		   delete[] p_mu[i];
		delete[] p_mu;
} // End parallel
		return act_A;
}




double pEnumeratorDouble::solveSVP(mat_ZZ& B, vec_ZZ& vec) {
	// Annealing mode was chosen
	if(Configurator::getInstance().do_annealing) {
		AnnealInfo<long double> ainfo;
		ainfo._number_of_random_bases = Configurator::getInstance().ann_bkz_instances;
		ainfo._number_of_annealing_threads = Configurator::getInstance().ann_annealing_threads;
		ainfo._number_of_parallel_reducing_threads = Configurator::getInstance().ann_parallel_reducing_threads;
		ainfo._number_of_different_bases = Configurator::getInstance().ann_num_different_bases;
		ainfo._time_per_node = Configurator::getInstance().ann_time_per_node; //4.145655e-8; //  3.082812277e-9

		SimulatedAnnealer<long double> annealer = SimulatedAnnealer<long double>(ainfo);
		annealer.anneal(B);
	}

	// Experimental use of fpllls strategizer functions
	if (Configurator::getInstance().ext_pruning_function_file.compare("NOT") != 0) {
		cout << "Reading pruning function from " << Configurator::getInstance().ext_pruning_function_file << endl;
		readPruning(Configurator::getInstance().ext_pruning_function_file);
	}

	// Several runs to test the difference
	double** p_mu = new double*[B.NumRows()];
	for(int i=0; i<B.NumRows(); i++)
		p_mu[i] = new double[B.NumRows()];

	double* p_bstar = new double[B.NumRows()];
	int* p_u = new int[B.NumRows()+2];
	double* p_prunfunc = new double[B.NumRows()+1];

	vec_ZZ randint;
	randint.SetLength(5);
	randint[0]=-2; randint[1]=-1; randint[2]=0; randint[3]=1; randint[4]=2;
	double act_A = -1;

	double tstart_oall = omp_get_wtime();

	mat_RR mu;
	vec_RR bstar;
	std::vector<mat_ZZ> Bcands;
	std::vector<int> Bcandsorder;

	int numt = 1;
	int procnum = 0;
	bool iterative_enumeration = Configurator::getInstance().iterative_enumeration;

	double BKZmaxtime = 9000;
	bool shortest_time_set = false;
        
	double target_a = -1;
	int dimen = B.NumCols();

	double theAval = Configurator::getInstance().Amax;
	if(theAval != numeric_limits<double>::max() && theAval > 0)	{
		target_a = theAval * theAval;
	}
	else if (Configurator::getInstance().do_gauss_estimate && Configurator::getInstance().gauss_estimate > 1e-5) {
		target_a = Configurator::getInstance().gauss_estimate;
	}
	else if(theAval < 0) {
	   target_a = calcGaussHeuristicChallenge<ZZ>(B, dimen);
	   target_a = target_a * target_a;
	}
	else
		target_a = p_bstar[0] * 1.00001;

	act_A = target_a;

	_conf = pruning_conf {Configurator::getInstance().enum_prune,
		Configurator::getInstance().prune_param};
	updatePruningFuncLoc(prunfunc.data(), _conf, act_A, dimen, 0, dimen-1);

#pragma omp parallel
#pragma omp single
		{
			numt = Configurator::getInstance().BKZinstances;
            if(numt<1) {
                numt = omp_get_num_threads();
                if(numt<1) numt=1;
            }
            
            else {
                numt= min(omp_get_num_threads(), Configurator::getInstance().BKZinstances);
            }
                

			if(numt<1) numt=1;

			cout << "Preparing preprocessing of " << numt << " threads." << endl;
			Bcands.reserve(numt);
			Bcands.resize(numt);
			Bcandsorder.reserve(numt);
			Bcandsorder.resize(numt);
		}
		srand (time(NULL));

        omp_set_num_threads(numt); // Use numt threads for all consecutive parallel regions
	    for(int trial=0; trial < Configurator::getInstance().trials; trial++)
	    {
			//double tprepstart = omp_get_wtime();

#pragma omp parallel num_threads(numt)
			{
			int tid = omp_get_thread_num();

			// Only numt thread should reduce in parallel
			if(tid < numt) {

			// Reset random basis
			Bcands[tid] = B;
			Bcandsorder[tid] = tid;
			if(omp_get_thread_num() == 0 && trial == 0 && Configurator::getInstance().test_original_basis) {
				// Keep basis unrandomized

			}
			else	{
				int seed = ((trial+1)*numt + (omp_get_thread_num()+1)*time(NULL));
				Bcands[tid] = randomizeMatrix(Bcands[tid], randint, Bcands[tid].NumCols(),seed);
			}

			// Do BKZ in parallel and generate a bunch
			//std::cout << "What reduction?" << std::endl;
			int prebeta=Configurator::getInstance().prebeta;
			if(Configurator::getInstance().dolll)
			{
				if(Configurator::getInstance().glob_beta <= 0 && prebeta <= 0)
				{
					double tstartlll = omp_get_wtime();
					std::cout << "LLL reducing basis" << std::endl;
					LLL_XD(Bcands[tid], Configurator::getInstance().glob_delta);
					double tendlll = omp_get_wtime();
					std::cout << "|b_0| of LLL reduced basis "
							<< lengthOfVector<ZZ, RR> (Bcands[tid][0])
							<< "(" << tendlll-tstartlll << "s)."
							<< std::endl;
				}

				else {
					bool old_auto = Configurator::getInstance().do_bkz_auto_tuning;
					for(int bval=Configurator::getInstance().glob_beta; bval<=Configurator::getInstance().glob_beta; bval++) {

						//DGB
						//std::cout << "Thread " << tid << " is BKZ-" << bval << " reducing basis" << std::endl;

						double tstartbkz = omp_get_wtime();
                        double tbkzexact_start = 1;
                        double tbkzexact_end = 0;
                        double tbkzheu_start = 1;
                        double tbkzheu_end = 0;

						// Alternative calculation in fpllls lib
						if(Configurator::getInstance().use_fplll_bkz)
						{
							// RAM based
							ZZ_mat<mpz_t> m_tored;
							ZZ_mat<mpz_t> U;

							m_tored.set_cols(B.NumCols());
							m_tored.set_rows(B.NumRows());

							for(int r=0; r < B.NumRows(); r++) {
								for(int c=0; c < B.NumCols(); c++) {
									mpz_t aa;
									mpz_init(aa);

									std::stringstream ssa;
									ssa << Bcands[tid][r][c];
									mpz_set_str( aa, ssa.str().c_str(),10);

									m_tored[r][c] = aa;
								}
							}

							vector<Strategy> strategies;


							// First run of pure BKZ with small beta
                            double bkzexacttime_start = omp_get_wtime();
							if(Configurator::getInstance().prebeta > 0)
							{
								//cout << "BKZ" << endl;
								BKZParam param20(prebeta, strategies);
	                            if(prebeta > 20) {
	                                //param20.flags |=  BKZ_MAX_TIME;// | BKZ_AUTO_ABORT; BKZ_VERBOSE;// |
	                                param20.flags |=  BKZ_AUTO_ABORT;
	                            }
	                            param20.delta = Configurator::getInstance().glob_delta;
	                            //param20.max_loops = 1;
								//param20.max_time = BKZmaxtime;
	                            tbkzexact_start = omp_get_wtime();
								bkz_reduction 	(&m_tored, &U,  param20, FT_DEFAULT, 0);
	                            tbkzexact_end = omp_get_wtime();
							}
                            double bkzexacttime_end = omp_get_wtime();
                            
                            if(!shortest_time_set)
                            {
                                shortest_time_set = true;
                                BKZmaxtime = 2.1 * (bkzexacttime_end - bkzexacttime_start);
                                cout << "Setting BKZ-Max-Time to " << BKZmaxtime << " s." << endl;
                            }

							// Second run with BKZ 2.0
                            tbkzheu_start = omp_get_wtime();
							if(Configurator::getInstance().glob_beta > 0)
							{
								strategies = load_strategies_json(strategy_full_path("default.json"));
								BKZParam paramval(bval, strategies);
								paramval.delta = Configurator::getInstance().glob_delta;
								paramval.flags |= BKZ_AUTO_ABORT;

	                            //for(int beta_i=4; beta_i<=Configurator::getInstance().glob_beta; beta_i+=4)
	                            {
	                            	//paramval.block_size = beta_i;
	                            	//bkz_reduction 	(&m_tored, &U,  paramval, FT_DEFAULT, 0);
	                            	bkz_reduction 	(&m_tored, &U,  paramval, FT_DEFAULT, 0);
	                            }
							}
							tbkzheu_end = omp_get_wtime();
                            
							for(int r=0; r < B.NumRows(); r++) {
								for(int c=0; c < B.NumCols(); c++) {
									std::stringstream ssb;
									ssb << m_tored[r][c];
									Bcands[tid][r][c] = conv<ZZ>(ssb.str().c_str());
								}
							}

						}

						else {
							Configurator::getInstance().do_bkz_auto_tuning = old_auto;
							LLL_XD(Bcands[tid], Configurator::getInstance().glob_delta);
							BKZ_FP(Bcands[tid], Configurator::getInstance().glob_delta, bval, 15);
						}

                        double tendbkz = omp_get_wtime();
						std::cout << "|b_0| of BKZ-" << prebeta << "/" << bval  << " reduced basis of thread "
                                << tid << " is: " 
								<< lengthOfVector<ZZ, RR> (Bcands[tid][0])
								<< "(" << tendbkz-tstartbkz << "s "
                                << "[" << tbkzexact_end - tbkzexact_start << " s / "
                                <<  tbkzheu_end - tbkzheu_start << "s ]"
                                << ")."                                
								<< std::endl;

						if(Configurator::getInstance().reducedfilepath!="" && numt==1)
						{
							std::ofstream basestream2;
							basestream2.open(Configurator::getInstance().reducedfilepath.c_str(), ios_base::out);
							//basestream2 << Bcands[tid];
							cout << "Exported the base" << endl;
						}
						else if(Configurator::getInstance().reducedfilepath!="" && numt>1) {
							cerr << "Reduced basis export not supported for multiple instances!" << endl;
						}
					}
				}
			}

			//writeBasisQuality<ZZ>(Bcands[tid], "");
			//exit(-1);

			//double gauss = calcGaussHeuristic<ZZ>(randbase);
			//	cout << "Gaussian heuristic for lambda_1(b): "
			//			<< sqrt(gauss) << endl;

			} // Only m threads should work
	        } // END PARALLEL


			// Sort by quality
			//sortBasesByQuality(Bcands, Bcandsorder);
			//writeBasisQuality<ZZ>(Bcands[0], "");

			// If wished for debug, the first basis is replaced by an external reduced one
			if(Configurator::getInstance().prereduced_base != "") {
				readRandomLattice(Bcands[0], Configurator::getInstance().prereduced_base);
			}

			//double tprepend = omp_get_wtime();
			//cout << "Time for preprocessing: " << tprepend-tprepstart << " seconds." << endl;

			// Now walk through all randomized bases and enumerate
			for(int tidi=0; tidi < (int)ceil(numt*1.0); tidi ++) {
				int tid = Bcandsorder[tidi];

				cout << "Processing base of thread " << tid << endl;
				procnum++;
				// Calculate Gram-Schmidt
				GSO(Bcands[tid], mu, bstar);

				double tstart1 = omp_get_wtime();

				/*** Begin former EnumTester ***/
				for(int i=0; i < B.NumCols(); i++) {
					p_u[i] = 0;
					p_bstar[i] = 0;
					p_prunfunc[i] = 0;
				}

				for(int i = 0; i < B.NumRows(); i++) {
					p_bstar[i] = NTL::conv<double>(bstar[i]);
					for(int j = 0; j < B.NumCols(); j++) {
						p_mu[i][j] = NTL::conv<double>(mu[i][j]);
					}
				}
				int dim = B.NumCols();

				act_A = EnumDouble(p_mu, p_bstar, p_u, 0, dim-1, dim, act_A, 0);

				NTL::Mat<long int> Bd; Bd.SetDims(B.NumRows(), B.NumCols());

				for(int i = 0; i < B.NumRows(); i++) {
					for(int j = 0; j < B.NumCols(); j++) {
						Bd[i][j] = NTL::conv<long int>(Bcands[tid][i][j]);
					}
				}

				NTL::Vec<double> sol;	sol.SetLength(B.NumCols());
				for(int ii=0; ii<B.NumCols(); ii++)
					sol[ii] = (double)p_u[ii];

				calcVectorLengthDouble<long int, double>(sol, Bd);
			    double tend1 = omp_get_wtime();
			    std::cout << "Runtime of standard enumeration: " <<
			    		tend1 - tstart1 << " s." << std::endl;

				/*** End former EnumTest ***/

			    // Standard behavior of Extreme Pruning
			    // End if a vector smaller than target found
			    if(act_A < target_a && !iterative_enumeration) {
			    	double tend_oall = omp_get_wtime();
			    	cout << "Finishing after " << trial+1 << " attemps in "
			    			<< tend_oall - tstart_oall << " s (processed "
							<< procnum << " bases)."
			    			<< endl;
			    	return act_A;
			    }

			    // Iterative behavior: Search as long as there is no better solution
			    // than minimum
			    else if(iterative_enumeration && act_A < target_a) {
			    	// Also the target value given by user was found

				    if(act_A < Configurator::getInstance().Aend * Configurator::getInstance().Aend) {
				    	double tend_oall = omp_get_wtime();
				    	cout << "Finishing iterative search after " << trial+1 << " attemps in "
				    			<< tend_oall - tstart_oall << " s (processed "
								<< procnum << " bases)."
				    			<< endl;
				    	return act_A;
				    }

				    // The search has to go on. A short vector was found, but not the target
				    trial = 0; // Reset the trial count
			    }

			} // End of iterating through all candidates of this round
	    }


	    double tend_oall = omp_get_wtime();
	    	cout << "Unsuccessful after " << Configurator::getInstance().trials << " attemps in "
	    			<< tend_oall - tstart_oall << " s (processed "
					<< procnum << " bases)."
	    			<< endl;
	    	return -1;

		delete[] p_bstar;
		delete[] p_u;
		delete[] p_prunfunc;
		for(int i=0; i<B.NumRows(); i++)
		   delete[] p_mu[i];
		delete[] p_mu;
		return act_A;

}

void pEnumeratorDouble::resetEnumerator() {
	candidates_queue->clear();
	//candidates_vec.clear();

	prunfunc.SetLength(_dim+2);
	_conf = pruning_conf {NO, 1.0001};

	cand_s = 0; // highest non-zero entry
	cand_t = 0; // actual visited level in tree

	for(int i=0; i<=_dim+1; i++) {
		cand_umin[i] = 0;
		cand_v[i] = 0;
		cand_l[i] = 0;
		cand_c[i] = 0;
		cand_Delta[i] = 0;
		cand_delta[i] = 1;

		prunfunc[i] = 0;
	}

	for(int i=0; i < _num_threads+2; i++) {
		for(int j=0; j <_dim+2; j++) {
			umin[i][j] = 0;
			v[i][j] = 0.0;
			l[i][j] = 0.0;
			c[i][j] = 0.0;
			Delta[i][j] = 0;
			delta[i][j] = 0;
		}
	}

	_cnt = 0;
	_cnt2 = 0;
}


double pEnumeratorDouble::EnumTuner(double** mu, double* bstar, int* u, int beta, int dim, double A,
		int vec_offset, int jj, int kk) {
	if(jj<0 || kk <0) {
		jj = (dim - beta) / 2;
		kk = jj + beta;
	}

	cout << "Tuning for dim: " << kk - jj + 1 << endl;

	// What is the target length of the short vector?
	if(Configurator::getInstance().do_gauss_estimate && Configurator::getInstance().gauss_estimate > 1e-5) {
		A = Configurator::getInstance().gauss_estimate;
	}

	else {
		A = A*1.0001;
	}

	_conf = pruning_conf {NO, 1.0001};
	updatePruningFuncLoc(prunfunc.data(), _conf, A, kk-jj + 1, jj, kk);

	//int min_height = 0.21 * beta;
	//int max_height = 0.30 * beta;

	double min_time = std::numeric_limits<double>::max();

	int max_threads = 1;

	int bestheight = -1;
	int bestthreads = -1;

#pragma omp parallel
#pragma omp single
	max_threads = omp_get_num_threads();

	int height = 1;
	for(int height_per = 170; height_per <= 250; height_per += 5) {

		if((int)round((double)(height_per * beta) / 1000) > (double)height)
			height = (int)round((double)(height_per * beta) / 1000);
		else
			continue;

		int act_thread = max_threads;
		//while(act_thread >= 1) {
			this->resetEnumerator();
			_conf = pruning_conf {NO, 1.0001};
			updatePruningFuncLoc(prunfunc.data(), _conf, A, kk-jj + 1, jj, kk);

			omp_set_num_threads(act_thread);
			for(int i=0; i<dim-2; i++)
				u[i]=0;


			double tstart = omp_get_wtime();

			pEnumeratorDouble::BurgerEnumerationDoubleParallelDriver(mu, bstar,
					u, jj, kk, dim, A, height, true);

			double tend = omp_get_wtime();

			//cout << "Height: " << height << " / threads: " << act_thread << " required "
			//		<< tend - tstart << " s." << endl;

			if(tend-tstart < min_time) {
				Configurator::getInstance().bkz_height = height;
				min_time = tend-tstart;
				bestheight = height;
				bestthreads = act_thread;
			}

			act_thread = floor(act_thread / 2);
		//}
	}

	 cout << "Best config is height: " << bestheight << " / threads: " << bestthreads << " required "
			<< min_time << " s." << endl;

	omp_set_num_threads(max_threads);
	Configurator::getInstance().do_bkz_auto_tuning = false;

	return 0;
}

double pEnumeratorDouble::EnumDouble(double** mu, double* bstar, int* u, int jj, int kk, int dim,
		double A, int vec_offset, bool is_bkz_enum) {

	cout << setprecision(5);
	int len = kk-jj + 1;

	double ret = 0.0;

	if(len > Configurator::getInstance().par_threshold) {
		int height = (int)round(0.22 * len);
		if(height < 1)
			height = 1;

		if(Configurator::getInstance().forced_height > 0 && !is_bkz_enum) {
			height = Configurator::getInstance().forced_height;
		}

		if(Configurator::getInstance().bkz_height > 0 && is_bkz_enum) {
			height = Configurator::getInstance().bkz_height;

		}

		ret = BurgerEnumerationDoubleParallelDriver(mu, bstar, u, jj, kk, dim+vec_offset, A, height, is_bkz_enum);

	}

	else {
		ret =  BurgerEnumerationDouble(mu, bstar, u, prunfunc, jj, kk, dim+vec_offset, A);
	}

	return ret;
}

double pEnumeratorDouble::BurgerEnumerationDoubleParallelDriver(double** mu, double* bstar,
		int* u, int jj, int kk, int dim, const double A, int candidate_height, bool is_bkz_enum) {

    _cand_cnt = 0;    
    _cnt2 = 0;
	VectorStorage::getInstance().cand_count = 0;
	MB::MBVec<int> startmin; startmin.SetLength(dim+2);
	MB::MBVec<int> ures; ures.SetLength(dim+2);
	for(int i=0; i <= dim+1; i++) {
		u[i] = 0;
		ures[i] = 0;
		startmin[i] = 0;
	}

	u[0] = 1;
	startmin[0] = 0;
	int serial_height = candidate_height;

	// Newest version
	bool candidates_left = true;
	bool below_thres = false;
	resetCandidateSearch();

	cand_t = kk-serial_height;
	//cand_t = 0;

	double Amin = prunfunc[jj] * 1.000;
	long long locnodeinit = 0;
	double enumtime_start = omp_get_wtime();

	cout << "Searching candidates in depth " << serial_height << endl;

	int ret = BurgerEnumerationCandidateSearch(mu, bstar,
	 		startmin, NULL, jj, kk-serial_height, kk, dim, locnodeinit, Amin);

	//int ret = BurgerEnumerationCandidateSearch(mu, bstar,
	//	 		startmin, NULL, jj, jj, serial_height, dim, locnodeinit, Amin);

	/*if(candidates_queue->checkForDuplicates())
		cout << "Duplicates!!!" << endl;
	else
		cout << "No duplicates." << endl;*/

	//DBG
	/*resetCandidateSearch();
	candidates_queue->clear();
	MB::MBVec<int> u_s; u_s.SetLength(dim+2);
	u_s[79] = 1;
	u_s[73] = -1;
	u_s[72] = -1;
	u_s[70] = -1;
	u_s[69] = 1;
	candidates_queue->push(u_s);
	candidates_queue->print(kk-serial_height);*/


	if(ret == 0) {
		candidates_left = false;
	}

	fflush(stdout);
	bool search_is_executed = false;
    bool early_out = false;
    int enumt = Configurator::getInstance().Enumthreads;
    omp_set_num_threads(enumt);
    cout << "Enumerating with " << Configurator::getInstance().Enumthreads << " threads." << endl; 
	long long locnodecnt = 0;

#pragma omp parallel shared(candidates_left) num_threads(enumt) reduction(+:locnodecnt)
	{
		bool do_candidate_search = false;
		bool do_enumeration = false;
		MB::MBVec<int> u_loc; u_loc.SetLength(dim+2);
		MB::MBVec<double> prunfuncloc; prunfuncloc.SetLength(dim+2);
		prunfuncloc = prunfunc;
		double Aret = Amin;

		while(true)
		{
			do_candidate_search=false;
			do_enumeration = false;

			// Check whether to execute the candidate search by this thread
#pragma omp single nowait
{
			below_thres = candidates_queue->isBelowThres();

			if(!search_is_executed) {
				if(candidates_left && below_thres) {
					do_candidate_search=true;
					search_is_executed=true;
				}
			}
}

			if(do_candidate_search) {
				int ret = BurgerEnumerationCandidateSearch(mu, bstar,
						startmin, NULL, jj, kk-serial_height, kk, dim, locnodecnt, Amin);
                        
				// Mark for all that there are no candidates and that no thread needs to do the search again
				if(ret == 0) {
					candidates_left = false;
				}
				else
					search_is_executed=false;
			}

			// Check whether there are still candidates to process
#pragma omp critical (QueueAccess)
{
			if(!candidates_queue->isEmpty()) {
				do_enumeration=true;
				candidates_queue->next(u_loc, kk-serial_height, kk);
                _cnt2++;
			}
}

			// Do enumeration if thread has something to do
			if(do_enumeration) {
				prunfuncloc = prunfunc;
				Aret = BurgerEnumerationDoubleRemainder(mu, bstar, u_loc,
						prunfuncloc.data(), jj, kk-serial_height, kk, dim, locnodecnt, Amin);

#pragma omp critical (ValuesUpdate)
				if(Aret < Amin)
				{
					Amin = Aret;
					ures = u_loc;
					if(!is_bkz_enum) {
					cout << "Found new minumun (length: " << (int)sqrt(Amin) << ")"
							<< ures << endl;
					}

					updatePruningFuncLoc(prunfunc.data(), _conf, Amin, kk-jj+1, jj, kk);
                    
                    if(Configurator::getInstance().enum_prune == EXTREME ||
                    		Configurator::getInstance().enum_prune == EVENLINEAR) {

                    	double enumtime_end = omp_get_wtime();
                    	cout << "Fastest solution: " << enumtime_end - enumtime_start << " s." << endl;
                        //cout << "Processed " << _cnt2 << " candidate-subtrees." << endl;
                    	//cout << "Speed:  " << (double)(locnodeinit + locnodecnt) / (enumtime_end - enumtime_start) << " nodes/s ("
                    	//		<< locnodeinit + locnodecnt << " nodes)." << endl;
                        early_out = true;
                    }
				}
			}


			if((!candidates_left && !do_enumeration) || early_out) {
#pragma omp cancel parallel
				break;
			}

		} // End main loop

	} // end parallel block

	for(int i=jj; i<=kk; i++) {
		u[i] = ures[i];
	}
	double enumtime_end = omp_get_wtime();

    cout << "Processed " << _cnt2 << " candidate-subtrees." << endl;
	cout << "Speed:  " << (double)(locnodeinit + locnodecnt) / (enumtime_end - enumtime_start) << " nodes/s ("
			<< locnodeinit + locnodecnt << " nodes)." << endl;
	return Amin;
}

void pEnumeratorDouble::resetCandidateSearch() {
	for(int i=0; i<=_dim+1; i++) {
		cand_umin[i] = 0;
		cand_v[i] = 0;
		cand_l[i] = 0;
		cand_c[i] = 0;
		cand_Delta[i] = 0;
		cand_delta[i] = 1;
	}

	candidates_queue->clear();
	//candidates_vec.clear();
}

/**

return: Integer indicating if candidates are left. If 0 is returned, no candidates are left in the tree
*/
int pEnumeratorDouble::BurgerEnumerationCandidateSearch(double** mu, double* bstarnorm,
		MB::MBVec<int>& u, const double* prun_func, const int min, int j, int k, int dim, 
		long long& locnodecnt,
		double Ain) {

	while (/*cand_t <= k*/true) {
			cand_l[cand_t] = cand_l[ cand_t + 1 ] +
					( (u[cand_t]) + cand_c[cand_t] ) * ( (u[cand_t]) + cand_c[cand_t] ) * bstarnorm[cand_t];

			if(cand_l[cand_t] < prunfunc[cand_t]) {
				locnodecnt++;
				if(cand_t > j) {
					cand_t = cand_t - 1;

					cand_r[cand_t] = std::max<int>(cand_r[cand_t], cand_r[cand_t+1]);

					for(int j = cand_r[cand_t+1]; j > cand_t + 1; j--) {
						cand_sigma[cand_t][j-1] = cand_sigma[cand_t][j] + u[j-1] * mu[j-1][cand_t];
					}
					cand_r[cand_t+1] = cand_t + 1;

					cand_c[cand_t] = cand_sigma[cand_t][cand_t+1];

					// Old: Recalculate everything
					//cand_c[cand_t] = muProdDouble(u.data(), mu, cand_t, k);

					u[cand_t] = cand_v[cand_t] = (optroundF(-cand_c[cand_t] - 0.5));

					/*if(abs(cand_sigma[cand_t][cand_t+1] - cand_c[cand_t]) < 1e-5) {
						cout << "Check" << endl;
					}

					else {
						cout << "cand_t: " << cand_t << endl;
						cout << cand_c[cand_t] << " vs. " << cand_sigma[cand_t][cand_t+1] << endl;
						exit(-1);
					}*/

					// For ZigZag
					cand_Delta[cand_t] = 0;
					if(u[cand_t] > -cand_c[cand_t]) {
						cand_delta[cand_t] = - 1;
					}
					else {
						cand_delta[cand_t] = 1;
					}
				}

				else {

#pragma omp critical (QueueAccess)
{
					candidates_queue->push(u, min, k); // Pushing vectors with full length
                    _cand_cnt++;
}

#ifdef VECDEBUG
					VectorStorage::getInstance().cand_count++;
#endif

					if(candidates_queue->isFull()) {
						return 1;
					}

					// Avoid visiting short vector twice
					// When bound is not updated

					cand_t = cand_t + 1;

					if(cand_t > k) {
						return 0;
					}

					cand_r[cand_t] = cand_t+1;

					cand_s = max(cand_s,cand_t);
					if (cand_t < cand_s) // If active: Half
						cand_Delta[cand_t] = -cand_Delta[cand_t];

					if(cand_Delta[cand_t] * cand_delta[cand_t] >= 0)
						cand_Delta[cand_t] += cand_delta[cand_t];

					u[cand_t] = cand_v[cand_t] + cand_Delta[cand_t];

				}
			}


			else {
				cand_t = cand_t + 1;

				if(cand_t > k) {
					return 0;
				}

				cand_r[cand_t] = cand_t+1;

				cand_s = max(cand_s,cand_t);

				if (cand_t < cand_s) // If active: Half
					cand_Delta[cand_t] = -cand_Delta[cand_t];

				if(cand_Delta[cand_t] * cand_delta[cand_t] >= 0)
					cand_Delta[cand_t] += cand_delta[cand_t];

				u[cand_t] = cand_v[cand_t] + cand_Delta[cand_t];
			}

		}

	return 0;
}


double pEnumeratorDouble::BurgerEnumerationDoubleRemainder(double** mu, double* bstarnorm,
		MB::MBVec<int>& u, double* prunefunc_in, int j, int k, int rel_len, int dim, long long& locnodecnt, double Ain) {

	const int myid = omp_get_thread_num();
	//double** sigma_ptr = sigma[myid];

	for(int i=0; i <_dim+2; i++) {
		for(int j=0; j < _dim+2; j++) {
			sigma[myid][i][j] = 0.0;
		}
		r[myid][i+1] = k+1;
	}

	MB::MBVec<double> sig_cache;
	sig_cache.SetLength(dim+3);

	int s;
	int t = j; // s required for zigzagpattern

	double A = 0;
	if(Ain == -1)
		A = bstarnorm[j] * 1.0001;
	else
		A = Ain;

	for(long ii=0; ii <= dim+1; ii++) {
		c[myid][ii] = l[myid][ii] = 0;
		v[myid][ii] = 0;
        Delta[myid][ii] = 0;
		delta[myid][ii] = 1;
		umin[myid][ii] = 0;
	}

	for(long ii=j; ii <= k; ii++) {
		u[ii] = 0.0;
	}

	// My start
	u[rel_len+1] = 0.0;

	int level = 0;
	long long xabssum = 0;
	double y;

	s = t = k;

	//NEWfor (level = rel_len; level >= j; level--) {
	for (level = rel_len; level >= k; level--) {
		if (l[myid][level+1] > prunefunc_in[level+1] && level+1 > k) {
			s = t = level + 1;
			return std::numeric_limits<double>::max(); // A shorter vector has been found by another thread
		}

		double ctmp = 0.0;
		ctmp = muProdDouble(u.data(), mu, level, rel_len);
		c[myid][level] = ctmp;
		y = u[level] + c[myid][level];
		l[myid][level] = l[myid][level+1] + y * y * bstarnorm[level];

		xabssum += (abs(u[level]));
	}

	for(int tloc=j; tloc<=k; tloc++) { // for all relevant values of t
		sig_cache[tloc] = 0;
		for(int vv=k ; vv < rel_len; vv++) { //min(k+1,tloc+1)
			sig_cache[tloc] += u[vv+1] * mu[vv+1][tloc];
		}
	}

	// Prepare the highest entry which may be changed
	u[t] = v[myid][t] = (optroundF(-c[myid][t] - 0.5));

	Delta[myid][t] = 0;

	if(u[t] > -c[myid][t]) {
		delta[myid][t] = - 1;
	}
	else {
		delta[myid][t] = 1;
	}

	// It is impossible that length is zero for non-zero vector
	// To prevent finding zero vector, then set [1, 0, 0, ... 0]
	//if (abs(l[myid][j]) < ) {
	if(xabssum < 1e-5) {
		u[j] = umin[myid][j] = 1.0;
		A = bstarnorm[j];
		l[myid][j] = prunefunc_in[j] /*A*/;
		s = t = j;
	}

	// Init the values of sigma and r
	for(int i=0; i <= k; i++) {
		sigma[myid][i][k+1] = sig_cache[i];
		int max_t = max<int>(s, i);

		for(int j=k+1; j > max_t+1; j--) {
			sigma[myid][i][j-1] = sigma[myid][i][j] + u[j-1]  * mu[j-1][i];
		}

		r[myid][i + 1] = max_t + 1;
	}


	//cout << u << "@t:" << t << endl;
	//int stop = 0;
	//cin >> stop;
	while(true) {

		//if(t > 67)
		//	cout << t << "\n";

		l[myid][t] = l[myid][ t + 1 ] + ( (u[t]) + c[myid][t] ) * ( (u[t]) + c[myid][t] ) * bstarnorm[t];

		if(l[myid][t] < prunefunc_in[t]) {
			locnodecnt++;
			if(t > j) {

				t = t - 1;
				//if(t==68)
				//	cout << u[68] << "@t:" << t << " (above)" << endl;

				r[myid][t] = std::max<int>(r[myid][t], r[myid][t+1]);

				for(int j = r[myid][t+1]; j > t + 1; j--) {
					sigma[myid][t][j-1] = sigma[myid][t][j] + u[j-1] * mu[j-1][t];
				}
				r[myid][t+1] = t + 1;

				c[myid][t] = sigma[myid][t][t+1];
				u[t] = v[myid][t] = (optroundF(-c[myid][t] - 0.5));

				// For ZigZag
				Delta[myid][t] = 0;

				if(u[t] > -c[myid][t]) {
					delta[myid][t] = - 1;
				}
				else {
					delta[myid][t] = 1;
				}
			}

			else {
				A = l[myid][t];
				updatePruningFuncLoc(prunefunc_in, _conf, l[myid][t], rel_len-j+1, j, k);

				for(long i=j; i <= rel_len; i++) {
						umin[myid][i] = u[i];
				}

				// Avoid visiting short vector twice
				// When bound is not updated
				t = t + 1;

				//if(t==68)
				//	cout << u[68] << "@t:" << t << " (below)" << endl;

				if(t > k)
					break;

				r[myid][t] = t+1;
				s = max(s,t);

				//if (t < s) If ative: half
					Delta[myid][t] = -Delta[myid][t];

				if(Delta[myid][t] * delta[myid][t] >= 0)
					Delta[myid][t] += delta[myid][t];

				u[t] = v[myid][t] + Delta[myid][t];
			}
		}


		else {
			t = t + 1;

			//if(t==68)
			//	cout << u[68] << "@t:" << t << " (below)" << endl;

			if(t > k)
				break;
			r[myid][t] = t+1;
			s = max(s,t);

			//if (t < s)
				Delta[myid][t] = -Delta[myid][t];

			if(Delta[myid][t] * delta[myid][t] >= 0)
				Delta[myid][t] += delta[myid][t];

			u[t] = v[myid][t] + Delta[myid][t];
			//if(t==68)
			//	cout << "u[68]: " << u[68] << endl;
		}

	}
	//cout << u << "@t (END):" << t << endl;
	for(long i=j; i <= rel_len; i++) {
		u[i] = umin[myid][i];
	}
	return A;
}


double pEnumeratorDouble::BurgerEnumerationDouble(double** mu, double* bstarnorm,
		int* u, MBVec<double> prunefunc_in, int j, int k, int dim, double Ain=-1) {

	int tid = omp_get_thread_num();
	long long nodecnt = 0;

	for(int i=0; i< _dim + 2; i++) {
		r[tid][i] = i;
		for(int j=0; j < _dim + 2; j++) {
			sigma[tid][i][j] = 0.0;
		}
	}

	int* umin = new int[dim + 2];
	int* v = new int[dim + 2];
	double* l = new double[dim + 2];
	double* c = new double[dim + 2];

	// To decide the zigzag pattern
	int* Delta = new int[dim + 2];
	int* delta = new int[dim + 2];

	int s;
	int t; // s d for zigzagpattern

	double A = bstarnorm[j] * 1.0001;
	if(Ain == -1)
		A = bstarnorm[j] * 1.0001;
	else
		A = Ain;

	for(long i=0; i <= dim + 1; i++) {
		c[i] = 0.0;
        l[i] = 0;
		Delta[i] = v[i] = 0.0;
		delta[i] = 1.0;
		umin[i] = u[i] = 0;
	}

	s = t = j;

	for(long i=j; i <= k+1; i++) {
		umin[i] = u[i] = 0;
	}
	u[j] = umin[j] = 1;


	double tstart = omp_get_wtime();
	while (true) {

		l[t] = l[ t + 1 ] + ( (u[t]) + c[t] ) * ( (u[t]) + c[t] ) * bstarnorm[t];

		if(l[t] < prunefunc_in[t]) {
			nodecnt++;
			if(t > j) {
				t = t - 1;

				r[tid][t] = std::max<int>(r[tid][t], r[tid][t+1]);

				for(int j = r[tid][t+1]; j > t + 1; j--) {
					sigma[tid][t][j-1] = sigma[tid][t][j] + u[j-1] * mu[j-1][t];
				}
				r[tid][t+1] = t + 1;

				c[t] = sigma[tid][t][t+1];

				// Activate to not use caching
				//c[t] = muProdDouble(u, mu, t, s/*dim-1*/);
                
                u[t] = v[t] = optroundF(-c[t] - 0.5);        

				// For ZigZag
				Delta[t] = 0;
				if(u[t] > -c[t]) {
					delta[t] = - 1;
				}
				else {
					delta[t] = 1;
				}
			}

			else {
				A = l[t];

				for(long i=j; i <= k; i++) {
					umin[i] = u[i];
				}

				// Avoid visiting short vector twice
				// When bound is not updated
				t = t + 1;

				if(t > k) {
					break;
				}

				r[tid][t] = t+1;

				s = max(s,t);

				if (t<s)
					Delta[t] = -Delta[t];

				if(Delta[t] * delta[t] >= 0) Delta[t] += delta[t];
					u[t] = v[t] + Delta[t];
			}
		}

		else {
			t = t + 1;


			if(t > k) {
				break;
			}

			r[tid][t] = t+1;

			s = max(s,t);

			if (t<s)
				Delta[t] = -Delta[t];

			if(Delta[t] * delta[t] >= 0) Delta[t] += delta[t];
			u[t] = v[t] + Delta[t];
		}

	}

	for(long i=j; i <= k; i++) {
		u[i] = umin[i];
	}
	double tend = omp_get_wtime();

	cout << "[Thread: " << tid << "]\t"
		<< "Time ENUM: " << (tend - tstart)
			<< " / Speed:  " << (double)(nodecnt) / (tend - tstart) << " nodes/s ("
			<< nodecnt << " nodes)." << endl;

	// clean memory
	delete[] umin;
	delete[] v;
	delete[] l;
	delete[] c;
	delete[] Delta;
	delete[] delta;

	return A;
}




inline double pEnumeratorDouble::muProdDouble (double* x, double** mu, const int t, const int s) {

	// Only small variant
	double res(0);
	for(int i = t + 1; i <= s; i++) {
		res += x[i] * mu[i][t];
//#pragma omp atomic
		//VectorStorage::getInstance().mult_cnt++;
	}
	return res;
}

inline double pEnumeratorDouble::muProdDouble (int* x, double** mu, const int t, const int s) {

	// Only small variant
	double res(0);
	for(int i = t + 1; i <= s; i++) {
		res += x[i] * mu[i][t];
//#pragma omp atomic
		//VectorStorage::getInstance().mult_cnt++;
	}
	return res;
}
