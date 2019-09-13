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
#include "boost/multi_array.hpp"

using namespace std;
using namespace NTL;
using namespace MB;

/**
    Assumes resulting value fits into short!
*/
inline int optroundF(const double &src) {
    return (32768 - (int)(32768. - src));
}

typedef boost::multi_array<double, 3> Dim3Array;

template<class T>
class MatrixDim3 {
	MatrixDim3 (int x, int y, int z) {
		_vals = new T[x*y*z];
		dim1=x;
		dim2=y;
		dim3=z;
	}

private:
	T* _vals;
	int dim1, dim2, dim3;
};

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
	double solveSVPMP(mat_ZZ& B, vec_ZZ& vec);

	double EnumDouble(double** mu, double* bstar, int* u, int jj, int kk, int dim, double A, int vec_offset=1, bool is_bkz_enum=false);
	double EnumTuner(double** mu, double* bstar, int* u, int beta, int dim, double A, int vec_offset=1, int jj=-1, int kk=-1);

	template <class FT>
	void testAnnealing(const mat_ZZ& B);

	long long get_cnt () {
		return _cnt;
	}
private:
	// Memberfunctions
	double BurgerEnumerationDoubleParallelDriver(double** mu, double* bstarnorm,
			int* u, int jj, int kk, int dim, const double A, int candidate_height=9, bool is_bkz_enum=false);

	void resetEnumerator();
	void resetCandidateSearch();

	int BurgerEnumerationCandidateSearch(double** mu, double* bstarnorm,
			MBVec<int>& u, const double* prun_func, const int min, int j, int k, int dim, long long& locnodecnt, double Ain=-1.0);

	double BurgerEnumerationDoubleRemainder(double** mu, double* bstarnorm,
			MBVec<int>& u, double* prunefunc_in, int j, int k, int rel_len, int dim, long long& locnodecnt, double Ain=-1.0 );

	double BurgerEnumerationDouble(double** mu, double* bstarnorm,
			int* u, MBVec<double> prunefunc_in, int j, int k, int dim, double Ain);

	inline double muProdDouble (double* x, double** mu, const int t, const int s);
    inline double muProdDouble (int* x, double** mu, const int t, const int s);

	// Membervariables for candidate search
	// General enumeration variables
	int* cand_umin;
	int* cand_v; // next coordinate
	double* cand_l; // actual costs
	double* cand_c; // centers of all levels in tree

	// To decide the zigzag pattern
	int* cand_Delta; // sign
	int* cand_delta; // offset
	int cand_s; // highest non-zero entry
	int cand_t; // actual visited level in tree


	// Membervariables for the remainder search
	// General enumeration variables
	int** umin;
	int** v;
	double** l;
	double** c;


	// To decide the zigzag pattern
	int** Delta;
	int** delta;

	// To execute the benchmark
	unsigned long long* node_cnt;
	double * tstart_bench;

	// To reduce number of multiplications
	int** r;
	//Dim3Array sigma;
	double*** sigma;

	MBVecQueue3<int>* candidates_queue;
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
