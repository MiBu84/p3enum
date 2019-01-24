/*
 * BurgerEnum.hpp
 *
 *  Created on: 25.05.2018
 *      Author: tuxed
 */

#ifndef SRC_BURGERENUM_HPP_
#define SRC_BURGERENUM_HPP_

#include <vector>
//#include "parallelDagdalen.hpp"

#include <iostream>
#include <math.h>
#include <omp.h>
#include <limits>
#include "NTL/mat_RR.h"
#include <immintrin.h>
#include "MBQueue.hpp"
#include <stdio.h>
#include <boost/progress.hpp>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>
#include <fstream>      // std::ofstream
#include <algorithm>    // std::reverse

#include "pEnumeratorDouble.hpp"
#include "Utils.hpp"

#ifdef USE_MB_TYPES
#define VECTYPE MB::Vec
#define MATTYPE MB::Mat
#include "MBMat.hpp"
using namespace MB;
#else
#define VECTYPE NTL::Vec
#define MATTYPE NTL::Mat
using namespace NTL;
#endif

using namespace std;

class pEnumeratorDouble;


MB::MBVecQueue<double> glob_candidates_queue(15);
int multcnt;

template<class ZT, class FT> void computeG_BKZConstant(long beta, long p, VECTYPE<FT>& constants) {
	   const double c_PI = 3.14159265358979323846264338328;
	   const double LogPI = 1.14472988584940017414342735135;

	   VECTYPE<FT> Log;
	   Log.SetLength(beta);

	   long i, j, k;
	   double x, y;

	   for (j = 1; j <= beta; j++)
	      Log(j) = log(double(j));

	   for (i = 1; i <= beta-1; i++) {
	      // First, we compute x = gamma(i/2)^{2/i}

		   k = i/2;

		      if ((i & 1) == 0) { // i even
		         x = 0;
		         for (j = 1; j <= k; j++)
		            x = x + Log(j);

		         x = x * (1/double(k));

		         x = exp(x);
		      }
		      else { // i odd
		         x = 0;
		         for (j = k + 2; j <= 2*k + 2; j++)
		            x = x + Log(j);

		         x = 0.5*LogPI + x - 2*(k+1)*Log(2);

		         x = x * (2.0/double(i));

		         x = exp(x);
		      }
			  // Second, we compute y = 2^{2*p/i}

			  y = -(2*p/double(i))*Log(2);
			  y = exp(y);
			  constants(i) = x*y/c_PI;
	   }
}

template<class ZT, class FT> void ComputeG_BKZThresh( const VECTYPE<FT>& bstar, long beta, VECTYPE<FT>& threshold,
		const VECTYPE<FT>& constants)
{
   double x = 0;

   for (long i = 1; i <= beta-1; i++) {
      x += log(bstar[i-1]);
      threshold(i) = exp(x/double(i))*constants(i);
      if (!IsFinite(&threshold(i))) threshold(i) = 0;
   }
}


template<class ZT, class FT> inline void updatePruningFuncLoc(VECTYPE<FT>& prunfunc, const FT A, const int dim) {
	int i = 0;
	for(auto it=prunfunc.begin(); it!=prunfunc.end();++it) {
		*it = A;
		//*it = A * min(1.0, 1.05*(dim - i) / (dim - 1)); // 1.05
		/*if(*it <= 1e-4) {
			cerr << "Pruning func too small: p[" << i << "] = "<< *it
					<< " for A=" << A << " / dim: "<< dim << endl;
		}*/
		i++;
	}
}




template<class ZT, class FT> FT muProdDouble (const VECTYPE<FT>& x, const MATTYPE<FT>& mu, const int t, const int s) {

	// Only small variant
	if(s - t < 512 || 1==1) {
		FT res(0);
		for(int i = t + 1; i <= s; i++) {
			res += x[i] * mu[t][i];
			//multcnt++;
		}
		return res;
	}

	int tp1 = t+1;
	int istart = tp1 - (tp1%4);
	int its = (s - istart) / 4;
	alignas(16) FT sum2[4];
	sum2[0]=0; sum2[1]=0; sum2[2]=0; sum2[3]=0;

#ifdef _USE_AVX2
	__m256d ymm_sum = _mm256_load_pd(sum2);

	for(int i = istart; i <= s-3; i+=4) {
		//const __m256d ymm_mu = _mm256_load_pd(mu.data(t) + i);
		const __m256d ymm_x = _mm256_load_pd(x.data() + i);
		//ymm_sum = _mm256_fmadd_pd (ymm_x, ymm_mu, ymm_sum);
		cout << "ERROR" << endl;
	}

	_mm256_store_pd(sum2,ymm_sum);

#else
	for(int i = istart; i <= s-3; i+=4) {
		sum2[0] += (x[i+0]) * mu[t][i+0];
		sum2[1] += (x[i+1]) * mu[t][i+1];
		sum2[2] += (x[i+2]) * mu[t][i+2];
		sum2[3] += (x[i+3]) * mu[t][i+3];
	}
#endif

	// Remainder
	for(int i = istart + its*4; i <= s; i++) {
		sum2[0] += (x[i]) * mu[t][i];
	}

	FT res3 = sum2[0]+sum2[1]+sum2[2]+sum2[3];
	return res3;
}

template<class ZT, class FT> int BurgerEnumerationCandidateSearch(const MATTYPE<ZT>& B, const MATTYPE<FT>& mu, const VECTYPE<FT>& bstarnorm,
		VECTYPE<FT>& u, const VECTYPE<FT>& prun_func, int j, int k, int dim, int mode=2) {
	// When subtrees are enumerated, then all possible sign-combinations must be checked
		static VECTYPE<FT> umin; umin.SetLength(dim + 1);//(k - j + 2);
		umin = u;
		static VECTYPE<FT> v; v.SetLength(dim + 1);//(k - j + 2);
		static VECTYPE<FT> l; l.SetLength(dim + 1);//(k - j + 2);
		static VECTYPE<FT> c; c.SetLength(dim + 1);//(k - j + 2);

		// To decide the zigzag pattern
		static VECTYPE<FT> Delta; Delta.SetLength(dim + 1);//(k - j + 2);
		static VECTYPE<FT> delta; delta.SetLength(dim + 1);//(k - j + 2);

		// For caching of muprods
		static MATTYPE<FT> sigma; sigma.SetDims(dim+1,dim+1);
		static VECTYPE<FT> r; r.SetLength(dim+1);

		static int s;
		static int t; // s required for zigzagpattern

		if(mode == 2)
		{
			for(long i=0; i <= B.NumCols(); i++) {
				c[i] = l[i] = 0.0;
				Delta[i] = v[i] = 0.0;
				delta[i] = 1.0;
			}
			// Set scalars
			s = t = j;

			// It is impossible that length is zero for non-zero vector
			// To prevent finding zero vector, then set [1, 0, 0, ... 0]
			/*if (l[0] == 0.0)
				u[j] = 1.0;*/

			umin = u;
		}

		while (t <= k) {
			l[t] = l[ t + 1 ] + ( (u[t]) + c[t] ) * ( (u[t]) + c[t] ) * bstarnorm[t];
			if(l[t] < prun_func[t]) {
				if(t > j) {
					t = t - 1;
					// New: Use cached values
					// Adapt all coefficients below
					/*r[t] = max(r[t], r[t+1]);
					for(int jj=r[t]; jj >= t+1; jj--) {
						sigma[jj][t] = sigma[jj+1][t] + u[jj] * mu[t][jj];
					}
					c[t] = sigma[t+1][t];*/

					// Old: Recalculate everything
					c[t] = muProdDouble<ZT, FT>(u, mu, t, dim-1);
					u[t] = v[t] = (ceil(-c[t] - 0.5));

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
					glob_candidates_queue.push(u);

					// Avoid visiting short vector twice
					// When bound is not updated
					t = t + 1;
					r[t] = t;

					s = max(s,t);

					//if(!fix_sign)
					//	Delta[t] = -Delta[t];
					if (t<s)
						Delta[t] = -Delta[t];

					if(Delta[t] * delta[t] >= 0) Delta[t] += delta[t];
					u[t] = v[t] + Delta[t];

					if(glob_candidates_queue.isFull()) {
						//cerr << "Queue is full." << endl;
						return 1;
					}
				}
			}


			else {
				t = t + 1;
				//r[t] = t;
				s = max(s,t);

				/*if(!fix_sign)
					Delta[t] = -Delta[t];*/
				if (t<s)
					Delta[t] = -Delta[t];

				if(Delta[t] * delta[t] >= 0) Delta[t] += delta[t];
				u[t] = v[t] + Delta[t];
			}

		}

		return 0;

}


template<class ZT, class FT> static FT BurgerEnumerationDouble(const MATTYPE<ZT>& B, const MATTYPE<FT>& mu, const VECTYPE<FT>& bstarnorm,
		VECTYPE<FT>& u, VECTYPE<FT>& prunefunc_in, int j, int k, int dim, int mode=0, FT Ain=-1) {

	// When subtrees are enumerated, then all possible sign-combinations must be checked
	bool fix_sign = false;
	if(j==0 && k==dim-1)
		fix_sign = true;

	VECTYPE<FT> umin; umin.SetLength(dim + 1);//(k - j + 2);
	VECTYPE<FT> v; v.SetLength(dim + 1);//(k - j + 2);
	VECTYPE<FT> l; l.SetLength(dim + 1);//(k - j + 2);
	VECTYPE<FT> c; c.SetLength(dim + 1);//(k - j + 2);

	// To decide the zigzag pattern
	VECTYPE<FT> Delta; Delta.SetLength(dim + 1);//(k - j + 2);
	VECTYPE<FT> delta; delta.SetLength(dim + 1);//(k - j + 2);

	// For caching of muprods
	MATTYPE<FT> sigma; sigma.SetDims(dim+1,dim+1);
	VECTYPE<ZT> r; r.SetLength(dim+1);

	int s;
	int t; // s required for zigzagpattern
	FT A = Ain;
	if(Ain == -1)
		A = std::numeric_limits<double>::max();

	r[0] = 0;
	for(long i=0; i <= B.NumCols(); i++) {
		c[i] = l[i] = 0.0;
		Delta[i] = v[i] = 0.0;
		delta[i] = 1.0;
		r[i] = i;

		for(long j=0; j <= B.NumCols(); j++) {
			sigma[i][j] = 0.0;
		}
	}

	s = t = j;

	// initial calculation of c and l if some parts of the tree
	// are already processed
	// fill the mu-prod-cache for values which are missing
	// fill the mu-prod-cache for values which are missing
	/*for(int tv=0; tv <= dim-1; tv++) {
		for(int ii=1; ii < dim-t; ii++) {
			sigma[tv][ii] = sigma[tv][ii-1] + u[dim-ii] * mu[tv][dim-ii];
		}
		c[tv] = sigma[tv][dim-t-1];
	}*/


	for(int i = j; i < dim; i++) {
		c[i] = muProdDouble<ZT, FT>(u, mu, i, dim-1);
	}

	for(int i=dim-1; i>=0; i--) {
		l[i] = l[ i + 1 ] + ( u[i] + c[i] ) * ( u[i] + c[i] ) * bstarnorm[i];
	}

	// It is impossible that length is zero for non-zero vector
	// To prevent finding zero vector, then set [1, 0, 0, ... 0]
	if (l[0] == 0.0 && j == 0.0) {
		cout << "Setting last entry to 1." << endl;
		u[0] = 1.0;
	}

	umin = u;

	// Set scalars
	multcnt=0;

	// Prepare SchnorrHoerner Pruning
	VECTYPE<FT> constants;
	VECTYPE<FT> threshold;
	constants.SetLength(k-j-1);
	threshold.SetLength(k-j-1);
	computeG_BKZConstant<ZT, FT>(k-j, 10, constants);
	ComputeG_BKZThresh<ZT, FT>(bstarnorm, k-j, threshold, constants);
	std::reverse(threshold.begin(), threshold.end());

	while (t <= k) {
		l[t] = l[ t + 1 ] + ( (u[t]) + c[t] ) * ( (u[t]) + c[t] ) * bstarnorm[t];
		FT eta = threshold[t];
		if(l[t] < prunefunc_in[t]/* - eta*/) {
			if(t > j) {
				t = t - 1;

				// New: Use cached values for [t]
				// Working backward version
				/*r[t] = max(r[t], r[t+1]);
				for(int ii=r[t]; ii>=t; ii--) {
					sigma[t][ii] = sigma[t][ii+1] + u[ii+1] * mu[t][ii+1];
					//multcnt++;
				}
				c[t] = sigma[t][t];
				r[t+1] = t;*/


				// Working forward version
				r[t] = max(r[t], r[t+1]);
				//for(int ii=1; ii < dim-t; ii++) {
				/*for(int ii=dim-r[t]-1; ii < dim-t; ii++) {
					sigma[t][ii] = sigma[t][ii-1] + u[dim-ii] * mu[t][dim-ii];
				}
				c[t] = sigma[t][dim-1-t];*/
				//r[t+1] = t;


				// Working in serial
				//r[t] = min(r[t-1], r[t]);
				//const FT* p_u = u.data() + r[t];
				//const FT* p_mu = mu[t].data() + r[t];
				//for(int ii=1; ii < dim-1; ii++) {
					//sigma[t][ii] = sigma[t][ii-1] + u[dim-1-ii] * mu[t][dim-1-ii];
					//sigma[t][ii] = sigma[t][ii-1] + (*p_u) * (*p_mu);
					//p_u--;
					//p_mu--;
					//multcnt++;
				//}
				//c[t] = sigma[t][dim-2];
				//r[t+1] = t;

				// Old: Recalculate everything
				c[t] = muProdDouble<ZT, FT>(u, mu, t, dim-1);

				u[t] = v[t] = (ceil(-c[t] - 0.5));

				// For ZigZag
				Delta[t] = 0;
				if((FT)(u[t]) > -c[t]) {
					delta[t] = - 1;
				}
				else {
					delta[t] = 1;
				}
			}

			else {
				if(mode <= 1) {
					updatePruningFuncLoc<ZT,FT>(prunefunc_in, l[t], dim);
					A = l[t];
					umin = u;
				}

				// Just find one short vector
				if(mode==1) {
					return A;
				}

				// Avoid visiting short vector twice
				// When bound is not updated
				t = t + 1;
				r[t] = t;


				s = max(s,t);

				if (!fix_sign)
					Delta[t] = -Delta[t];
				else if (fix_sign && t<s)
					Delta[t] = -Delta[t];

				if(Delta[t] * delta[t] >= 0) Delta[t] += delta[t];

				u[t] = v[t] + Delta[t];
			}
		}


		else {
			t = t + 1;
			r[t] = t;
			s = max(s,t);

			if (!fix_sign)
				Delta[t] = -Delta[t];
			else if (fix_sign && t<s)
				Delta[t] = -Delta[t];

			if(Delta[t] * delta[t] >= 0) Delta[t] += delta[t];

			u[t] = v[t] + Delta[t];
		}

	}

	u = umin;
	return A;

}

template<class ZT, class FT> int BurgerEnumerationDoubleStandardDriver(const NTL::Mat<NTL::ZZ>& B, const NTL::Mat<NTL::RR>& mu, const NTL::Vec<NTL::RR>& bstar) {

	MATTYPE<ZT> Bd; Bd.SetDims(B.NumRows(), B.NumCols());
	MATTYPE<FT> mud; mud.SetDims(mu.NumRows(), mu.NumCols());
	VECTYPE<FT> bstard; bstard.SetLength(bstar.length());

	for(int i = 0; i < B.NumRows(); i++) {
		bstard[i] = NTL::conv<FT>(bstar[i]);
		for(int j = 0; j < B.NumCols(); j++) {
			Bd[i][j] = NTL::conv<long int>(B[i][j]);
			mud[j][i] = NTL::conv<FT>(mu[i][j]);
		}
	}

	VECTYPE<FT> u; u.SetLength(B.NumCols()+1);//(k - j + 2);
	for(int i=0; i <= B.NumCols(); i++) {
		u[i] = 0;
	}

	double tstart1 = omp_get_wtime();

	VECTYPE<FT> prunfunc;	prunfunc.SetLength(B.NumCols());
	updatePruningFuncLoc<ZT, FT>(prunfunc, bstard[0], B.NumCols());

	BurgerEnumerationDouble<ZT, FT>(Bd, mud, bstard, u , prunfunc, 0, B.NumCols()-1, B.NumCols(), 0, bstard[0]);
	calcVectorLengthDouble<ZT, FT>(u, Bd);
    double tend1 = omp_get_wtime();
    std::cout << "Runtime of standard enumeration: " <<
    		tend1 - tstart1 << " s." << std::endl;

	return 0;
}

template<class ZT, class FT> int BurgerEnumerationDoubleParallelDriver(const NTL::Mat<NTL::ZZ>& B, const NTL::Mat<NTL::RR>& mu, const NTL::Vec<NTL::RR>& bstar
		, int candidate_height=2) {
	int dim = B.NumCols();
	MATTYPE<ZT> Bd; Bd.SetDims(B.NumRows(), B.NumCols());
	MATTYPE<FT> mud; mud.SetDims(mu.NumRows(), mu.NumCols());
	VECTYPE<FT> bstard; bstard.SetLength(bstar.length());

	for(int i = 0; i < B.NumRows(); i++) {
		bstard[i] = NTL::conv<FT>(bstar[i]);
		for(int j = 0; j < B.NumCols(); j++) {
			Bd[i][j] = NTL::conv<ZT>(B[i][j]);
			mud[j][i] = NTL::conv<FT>(mu[i][j]);
		}
	}

	VECTYPE<FT> u; u.SetLength(B.NumCols()+1);//(k - j + 2);
	for(int i=0; i <= B.NumCols(); i++) {
		u[i] = 0;
	}

#ifdef FILEOUTPUT
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y--%H-%M-%S");
    string filename = oss.str();
    filename.append("Dim");
    filename.append(to_string(B.NumCols()));
    filename.append(".txt");
    cout << "Protocol-File: " << filename << std::endl;

    std::ofstream ofs;
    ofs.open (filename, std::ofstream::out);
    ofs << "Starting with A=" << bstard[0] << endl;
    ofs.close();
#endif

	int serial_height = candidate_height;
	cout << "Old A: " << bstard[0] << endl;
	//FT Amin = BurgerEnumerationDouble<ZT, FT>(Bd, mud, bstard, u , bstard[0], 0, B.NumCols()-1, B.NumCols(), 1);
	VECTYPE<FT> ures = u;
	FT Amin = bstard[0]*1.05;

	VECTYPE<FT> prunfunc;	prunfunc.SetLength(B.NumCols());
	updatePruningFuncLoc<ZT, FT>(prunfunc, Amin, B.NumCols());
	cout << "New A: " << Amin << " for u " << u << endl;

	for(int ii=0; ii <= B.NumCols()-serial_height; ii++) {
		u[ii] = 0;
	}

	double tstart = omp_get_wtime();

	// Newest version
	volatile bool candidates_left = true;
	bool below_thres = false;
	int mode = 2;

	// Do initial Candidate search
	printf("Starting initial candidate search.\n");

	// Testvec which is the shortest
	VECTYPE<FT> testvec;
	testvec.SetLength(B.NumCols() + 1);

	VECTYPE<FT> startmin;
	startmin.SetLength(B.NumCols() + 1);

	for(int ii=0; ii <= B.NumCols(); ii++) {
		startmin[ii] = 0;
		testvec[ii] = 0;
	}

	//testvec[39] = 0;
	//testvec[38] = 0;
	//testvec[36] = 1;
	//testvec[35] = -1;

	int ret = BurgerEnumerationCandidateSearch<ZT, FT>(Bd, mud, bstard,
			startmin , prunfunc, B.NumCols()-serial_height-1, B.NumCols()-1, B.NumCols(), mode);

	//glob_candidates_queue.push(testvec);
	//int ret = 0;

	printf("Finishing initial candidate search.\n");
	printf("%d candidates found.\n", glob_candidates_queue.elemsSize());

	mode++;
	if(ret == 0) {
		candidates_left = false;
		printf("All candidates processed initially.\n");
	}
	fflush(stdout);

	bool search_is_executed = false;

#pragma omp parallel shared(candidates_left)
	{
		bool do_candidate_search = false;
		bool do_enumeration = false;
		VECTYPE<FT> u_loc;
		VECTYPE<FT> prunfuncloc; prunfuncloc.SetLength(B.NumCols());
		prunfuncloc = prunfunc;
		FT Aret = Amin;

		while(true)
		{
			do_candidate_search=false;
			do_enumeration = false;

			// Check whether to execute the candidate search by this thread
#pragma omp single nowait
{
			below_thres = glob_candidates_queue.isBelowThres();
			if(!search_is_executed) {
				if(candidates_left && below_thres) {
					do_candidate_search=true;
					search_is_executed=true;
				}
			}
}

			if(do_candidate_search) {
				int ret = BurgerEnumerationCandidateSearch<ZT, FT>(Bd, mud, bstard,
											u, prunfunc, B.NumCols()-serial_height-1, B.NumCols()-1, B.NumCols(), mode);
				//printf("Thread %d releases the search.\n", omp_get_thread_num());
				// Mark for all that there are no candidates and that no thread needs to do the search again
				if(ret == 0) {
					candidates_left = false;
					printf("All candidates processed by thread %d.\n", omp_get_thread_num());
					fflush(stdout); // Prints to screen or whatever your standard out is
				}
				else
					search_is_executed=false;
			}

			// Check whether there are still candidates to process
#pragma omp critical (QueueAccess)
{
			if(!glob_candidates_queue.isEmpty()) {
				do_enumeration=true;
				u_loc = glob_candidates_queue.next();
				u_loc[B.NumCols()-serial_height-1] = 0;
			}
}

			// Do enumeration if thread has something to do
			if(do_enumeration) {
				//cout << "Max t = " << B.NumCols()-serial_height-1 << endl;
				prunfuncloc = prunfunc;
				Aret = BurgerEnumerationDouble<ZT, FT>(Bd, mud, bstard, u_loc,
						prunfuncloc, 0, B.NumCols()-serial_height-1, B.NumCols(), 0, Amin);
#pragma omp critical (ValuesUpdate)

				if(Aret < Amin)
				{
					Amin = Aret;
					ures = u_loc;
					updatePruningFuncLoc<ZT,FT>(prunfunc, Amin, dim);
					cout << "New A=" << Amin << " / coeffs: " << u_loc << endl;

#ifdef FILEOUTPUT
				    ofs.open (filename, std::ofstream::out | std::ofstream::app);
				    ofs << "New A=" << Amin << " / coeffs: " << u_loc << endl;
				    ofs.close();
#endif
				}
			}


			if(!candidates_left && !do_enumeration)
				break;

		} // End main loop

	} // end parallel block


#ifdef FILEOUTPUT
				    ofs.open (filename, std::ofstream::out | std::ofstream::app);
				    ofs << "Finished successfully." << endl;
				    ofs.close();
#endif

	cout << Amin << endl;

	calcVectorLengthDouble<ZT,FT>(ures, Bd);


    double tend = omp_get_wtime();
    std::cout << "Runtime of enumeration: " <<
    		tend - tstart << " s." << std::endl;


	return 0;

	//auto rng = std::default_random_engine {};
	//std::shuffle(std::begin(glob_candidates), std::end(glob_candidates), rng);


//
	//for(size_t i=0; i<glob_candidates.size(); i++) {


		//Aret = BurgerEnumerationDouble<ZT, FT>(Bd, mud, bstard, glob_candidates[i] , Amin, 0, B.NumCols()-serial_height - 1, B.NumCols(), 0);

		//cout << glob_candidates[i] << " / " << glob_candidates_queue.next() << endl;
		//cout <<  glob_candidates_queue.next() << endl;


	//}

}

template<class ZT, class FT> double BurgerEnumerationDoubleTester(const NTL::Mat<NTL::ZZ>& B, const NTL::Mat<NTL::RR>& mu,
		const NTL::Vec<NTL::RR>& bstar, double minA=-1 ) {
	double** p_mu = new double*[B.NumRows()];
	for(int i=0; i<B.NumRows(); i++)
		p_mu[i] = new double[B.NumRows()];

	double* p_bstar = new double[B.NumRows()];
	double* p_u = new double[B.NumRows()+2];
	double* p_prunfunc = new double[B.NumRows()+1];

	for(int i=0; i < B.NumCols(); i++) {
		p_u[i] = 0;
		p_bstar[i] = 0;
		p_prunfunc[i] = 0;
	}

	for(int i = 0; i < B.NumRows(); i++) {
		p_bstar[i] = NTL::conv<FT>(bstar[i]);
		for(int j = 0; j < B.NumCols(); j++) {
			p_mu[i][j] = NTL::conv<FT>(mu[i][j]);
		}
	}

	double tstart1 = omp_get_wtime();

	int dim = B.NumCols();
	pEnumeratorDouble pener = pEnumeratorDouble(dim);

	double inA = 0;
	if(minA > 0)
		inA = minA;
	else
		inA = p_bstar[0] * 1.001;


	double thelen = pener.EnumDouble(p_mu, p_bstar, p_u, 0, dim-1, dim, inA, 0);

	MATTYPE<ZT> Bd; Bd.SetDims(B.NumRows(), B.NumCols());

	for(int i = 0; i < B.NumRows(); i++) {
		for(int j = 0; j < B.NumCols(); j++) {
			Bd[i][j] = NTL::conv<long int>(B[i][j]);
		}
	}

	VECTYPE<FT> sol;	sol.SetLength(B.NumCols());
	for(int ii=0; ii<B.NumCols(); ii++)
		sol[ii] = p_u[ii];


	calcVectorLengthDouble<ZT, FT>(sol, Bd);
    double tend1 = omp_get_wtime();
    std::cout << "Runtime of standard enumeration: " <<
    		tend1 - tstart1 << " s." << std::endl;

	delete[] p_bstar;
	delete[] p_u;
	delete[] p_prunfunc;

	for(int i=0; i<B.NumRows(); i++)
	   delete[] p_mu[i];
	delete[] p_mu;
	return thelen;

}

#endif /* SRC_BURGERENUM_HPP_ */
