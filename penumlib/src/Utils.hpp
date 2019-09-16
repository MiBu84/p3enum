/*
 * Utils.hpp
 *
 *  Created on: 27.07.2018
 *      Author: mburger
 */

#ifndef SRC_UTILS_HPP_
#define SRC_UTILS_HPP_



#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include <sys/mman.h>
#include <unistd.h>
#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <tgmath.h>

//#include <boost/math/special_functions/gamma.hpp>

#include <NTL/RR.h>
#include <NTL/LLL.h>

#include "pruningfunc.hpp"
#include "Configurator.hpp"
#include "Utils.hpp"

using namespace std;
using namespace NTL;

NTL_CLIENT

Pruning to_prune(const std::string& s);
int readConfig(string path);
int readPruning(string path);
int readRandomLattice(NTL::mat_ZZ& B, std::string filepath);
int writeNTLLatticeToFile(const NTL::mat_ZZ& B, std::string filepath);

mat_RR GSO(const mat_ZZ& input, mat_RR& mu, vec_RR& c);
mat_ZZ randomizeMatrix(const mat_ZZ A, const vec_ZZ randoms, const int dim, const int seed);
void searchShortestSeedCandidates(int dim, int startseed, int endseed);
void sortBasesByQuality(const std::vector<mat_ZZ>& Bcands, std::vector<int>& order);

template<class T> void subMatrix(const mat_ZZ& in_mat, mat_ZZ& out_mat, int min, int max) {
	long int cols = max - min + 1;
	out_mat.SetDims(in_mat.NumRows(), cols);

	int i = 0;
	//int j = 0;

	for(int r = 0; r < in_mat.NumRows(); r++) {
		i = 0;
		for(int c=min; c < max; c++) {
			out_mat[r][i] = in_mat[r][c];
			i++;
		}
	}
}


template<class T>
string convert2Str(T &i) {
        string s;
        std::stringstream ss(s);
        ss << i;
        return ss.str();
        }

template <class T>
void convertFromStr(T &value, const std::string &s) {
        std::stringstream ss(s);
        ss >> value;
        }


template <class T>
double calcGaussHeuristicChallenge (const NTL::Mat<T>& _inputbasis, int n) {
	NTL::Mat<T> basis_trans; basis_trans.SetDims(_inputbasis.NumCols(), _inputbasis.NumRows());
	transpose(basis_trans,_inputbasis);

	//T d = abs(determinant(_inputbasis));
	T d = SqrRoot(determinant(basis_trans * _inputbasis));

	double det = conv<double>(d);
	det = pow(det, 1.0/double(n));

	double arg = ((double)n/2.0 + 1.0);
	double gam = exp(lgamma(arg));
	double frac = pow(gam, 1/double(n)) / sqrt(M_PI);

	double fpA = 1.05 * frac * det;
	//cout << "New Gauss is: "<< fpA << endl;

	/*double zahler = pow(exp(lgamma(arg)), 1.0/(double)n);
	double nenner = sqrt(M_PI);
	double volS = zahler/nenner;

	double det2 = conv<double>(d);
	det2 = pow(det2, 1.0/double(n));

	fpA = 1.05 * volS * det2;*/

	cout << "NewNew Gauss is: "<< fpA << endl;

	return fpA;
}

template <class T>
double estimateShortBase(const mat_ZZ& B, int start, int end) {
	mat_ZZ shortbase;
	subMatrix<double>(B, shortbase, start, end);
	double fpA = calcGaussHeuristicChallenge(shortbase, end-start+1);
	return fpA;
}

template <class T>
double calcGaussHeuristic (const NTL::Mat<T>& _inputbasis ) {
	int dim = _inputbasis.NumCols();

	NTL::Mat<T> basis = _inputbasis;
	T d = determinant(basis);
	RR dr;
	d = SqrRoot(abs(d));
	convertFromStr(dr, convert2Str(d)) ;
	RR dexp;
	double dp=double(2)/double(dim);
	convertFromStr(dexp, convert2Str(dp));
	dr = pow(dr, dexp);
	double dd;
	convertFromStr(dd, convert2Str(dr));
	double estSVL = pow(exp(gamma(dim/2+1)), 1/double(dim)) / sqrt(M_PI) * dd;
	double gauss = estSVL*estSVL;

	double fpA = gauss * 1.05;
	return fpA;
}



// Calculate length
template<class ZT, class FT> FT calcVectorLengthDouble(const NTL::Vec<FT>& coeffs, const NTL::Mat<ZT>& B) {
	//cout << coeffs << endl;
	//cout << B << endl;

	int d = B.NumCols();
	cout << "coefficients: [ ";
	for(int i=0; i < d; i++) {
		cout << coeffs[i] << " ";
	}
	cout << "]" << endl;

	FT len; len = 0;
	NTL::Vec<ZT> shortvec; shortvec.SetLength(d);

	for(int c = 0; c < d; c++) {
		shortvec[c] = 0;
	}

	std::stringstream vec_ss;
	for(int c = 0; c < d; c++) {
		for(int r = 0; r < d; r++) {
			shortvec[c] += NTL::conv<ZT>(coeffs[r]) * (B[r][c]);
		}
		vec_ss << shortvec[c] << " ";

	}

	for(int i=0; i < d; i++) {
		len += NTL::conv<FT>(shortvec[i] * shortvec[i]);
	}

	cout << "Shortest Vector (Len: "
			<< sqrt(len)
			<< "): ["
			<< vec_ss.str()
			<< "]" << endl;


	//cout << "Len: " << sqrt(len) << endl;
	cout << endl;
	return len;
}

template <class T> void printVector (T* vec, int len) {
	for(int i=0; i < len; i++) {
		cout << "[" << i << "]" << vec[i] << "   ";
	}
	cout << endl;
}

template <typename VT, typename FT> FT lengthOfVector(NTL::Vec<VT> vec) {
	FT len;
	len = 0.0;

	for(auto it = vec.begin(); it != vec.end(); it++) {
		len += conv<FT> (*it) * conv<FT> (*it);
	}
	len = NTL::sqrt(len);
	return len;
}



template <typename ZT> void writeBasisQuality(const NTL::Mat<ZT>& basis, std::string file) {
	mat_RR mu;
	vec_RR bstar;
	GSO(basis, mu, bstar);


	for(int i=0; i<basis.NumRows(); i++) {
		RR len = bstar[i];
		cout << i << "\t" << len << endl;
	}
}



#endif /* SRC_UTILS_HPP_ */
