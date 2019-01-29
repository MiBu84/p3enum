/*
 * Enumerator.hpp
 *
 *  Created on: 22.06.2018
 *      Author: mburger
 */

#ifndef SRC_PENUMERATORDOUBLE_HPP_
#define SRC_PENUMERATORDOUBLE_HPP_

#include <NTL/LLL.h>

#include "MBQueue.hpp"
#include "MBVec.hpp"
#include "pruningfunc.hpp"

using namespace std;
using namespace NTL;
using namespace MB;

class pEnumeratorDouble {
public:
	pEnumeratorDouble();
	pEnumeratorDouble(int dim);

	~pEnumeratorDouble();

	pEnumeratorDouble( const pEnumeratorDouble& other ) {
		cout << "pEnumeratorDouble( const pEnumeratorDouble& other )" << endl;
	}
	pEnumeratorDouble( pEnumeratorDouble& other ) {
		cout << "pEnumeratorDouble( const pEnumeratorDouble& other )" << endl;
	}

	double solveSVP(mat_ZZ& B, vec_ZZ& vec);
	double EnumDouble(double** mu, double* bstar, double* u, int jj, int kk, int dim, double A, int vec_offset=1, bool is_bkz_enum=false);
	double EnumTuner(double** mu, double* bstar, double* u, int beta, int dim, double A, int vec_offset=1, int jj=-1, int kk=-1);

	long long get_cnt () {
		return _cnt;
	}
private:
	// Memberfunctions
	double BurgerEnumerationDoubleParallelDriver(double** mu, double* bstarnorm,
			double* u, int jj, int kk, int dim, const double A, int candidate_height=9, bool is_bkz_enum=false);

	void resetEnumerator();
	void resetCandidateSearch();

	double BurgerEnumerationDoubleRemainder(double** mu, double* bstarnorm,
			MBVec<double>& u, double* prunefunc_in, int j, int k, int rel_len, int dim, long long& locnodecnt, double Ain=-1.0 );

	double BurgerEnumerationDouble(double** mu, double* bstarnorm,
			double* u, MBVec<double> prunefunc_in, int j, int k, int dim, double Ain);

	inline double muProdDouble (double* x, double** mu, const int t, const int s);

	int BurgerEnumerationCandidateSearch(double** mu, double* bstarnorm,
			MBVec<double>& u, const double* prun_func, const int min, int j, int k, int dim, long long& locnodecnt, double Ain=-1.0);



	// Membervariables for candidate search
	// General enumeration variables
	double* cand_umin;
	double* cand_v; // next coordinate
	double* cand_l; // actual costs
	double* cand_c; // centers of all levels in tree

	// To decide the zigzag pattern
	double* cand_Delta; // sign
	double* cand_delta; // offset
	int cand_s; // highest non-zero entry
	int cand_t; // actual visited level in tree


	// Membervariables for the remainder search
	// General enumeration variables
	double** umin;
	double** v;
	double** l;
	double** c;


	// To decide the zigzag pattern
	double** Delta;
	double** delta;

	// To execute the benchmark
	unsigned long long* node_cnt;
	double * tstart_bench;

	// To reduce number of multiplications
	//double** r;
	//volatile double*** sigma;

	MBVecQueue3* candidates_queue;
	//std::vector<MBVec<double> > candidates_vec;

	// General membervars
	int _dim;
	long long _cnt;
	long long _cnt2;
    long long _cand_cnt;
	int _num_threads;

	MBVec<double> prunfunc; // Array with pruning function
	pruning_conf _conf;
};


#endif /* SRC_PENUMERATORDOUBLE_HPP_ */
