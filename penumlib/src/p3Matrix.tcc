/*
 * p3Matrix.cpp
 *
 *  Created on: 08.06.2020
 *      Author: mburger
 */


#include "p3Matrix.hpp"
#include "Utils.hpp"

#include <iostream>
#include <algorithm>
#include <limits>
#include <cassert>
#include <string>
#include <cstring>

namespace p3 {
template<class T, class FT>
	p3Matrix<T, FT>::p3Matrix() {

	}

template<class T, class FT>
	p3Matrix<T, FT>::p3Matrix(int r, int c) {

	}

template<class T, class FT>
	p3Matrix<T, FT>::p3Matrix (const NTL::Mat<T>& mat) {
		this->_mat = mat;
		
	}

template<class T, class FT>
	void p3Matrix<T, FT>::BKZ () {	
}

template<class T, class FT>
	void p3Matrix<T, FT>::LLL_FP (double delta_) {
	cout << _mat << endl;
	//LLL_RR(_mat, delta_);
	//cout << _mat << endl;
	//return;

	FT delta = conv <FT>(delta_);
	_fmat = conv <Mat<FT>> (_mat);
	GSO();
	int k=1;

	while (k < getDim()) {
		sizeReduce(k, _fmat);

		if (_bstarabs[k] >= (delta - _mu[k][k-1] * _mu[k][k-1]) * _bstarabs[k-1] ) {
			k = k + 1;
		}

		else {
			swapBaseVectors(k, k-1);
			GSO();
			k = std::max(1,k-1);
		}
	}
}

template<class T, class FT>
	void p3Matrix<T, FT>::S2LLL (double eta) {	
	cout << _mat << endl;
	
	_fmat = conv <Mat<FT>> (_mat);
	GSO();
	
	int k = 1;
	while (k < getDim() ) {
		sizeReduce(k, _fmat);
		//std::cout << k << std::endl;

		// Compute i <-argmax
		int i = -1;
		std::vector<FT> S_jk;//NTL::conv<FT>(std::numeric_limits<double>min());
		argmax_S_lk(_fmat, k, S_jk, i);

		std::cout << S_jk[i] << std::endl;
		std::cout << calc_SS_B() << std::endl;

		if (S_jk[i] <= (1-eta) * calc_SS_B()) {
			k++;
		}

		else {
			sigma_i_k(i,k,_fmat);
			GSO();
			k = std::max(i, 1);
		}
	}
	cout << _mat << endl;
}

template<class T, class FT>
	void p3Matrix<T, FT>::sizeReduce(int k, Mat<FT>& fmat) {	
	
	for(int j =  k-1; j >= 0; j--) {
		FT q_j = round( _mu[k][j] );

		std::string strval = convert2Str<FT>(q_j);
		T q_j_i; convertFromStr<T>(q_j_i, strval);

		//fmat[k] = fmat[k] - q_j * fmat[j];
		_mat[k] = _mat[k] - q_j_i * _mat[j];
		_fmat = conv <Mat<FT>> (_mat);
		GSO();
	}
}

template<class T, class FT>
void p3Matrix<T, FT>::swapBaseVectors (int idx1, int idx2) {
	Vec<FT> vec1 = _fmat[idx1];
	Vec<FT> vec2 = _fmat[idx2];
	_fmat[idx1] = vec2;
	_fmat[idx2] = vec1;

	Vec<T> vec1i = _mat[idx1];
	Vec<T> vec2i = _mat[idx2];
	_mat[idx1] = vec2i;
	_mat[idx2] = vec1i;

	return;
}

template<class T, class FT>
int p3Matrix<T, FT>::sigma_i_k(int i, int k, Mat<FT>& fmat) {
	Vec<FT> row_back = fmat[i];
	Vec<FT> row_back2 = fmat[i];

	Vec<FT> row_k = fmat[k];

	fmat[i] = std::move(row_k);

	for(int pos=i; pos < k; pos++) {
		row_back2 = fmat[pos+1];
		fmat[pos+1] = row_back;
		row_back = row_back2;
	}


	Vec<T> row_backi = _mat[i];
	Vec<T> row_back2i = _mat[i];

	Vec<T> row_ki = _mat[k];

	_mat[i] = std::move(row_ki);

	for(int pos=i; pos < k; pos++) {
		row_back2i = _mat[pos+1];
		_mat[pos+1] = row_backi;
		row_backi = row_back2i;
	}

	return 0;
}

template<class T, class FT>
	void p3Matrix<T, FT>::GSO() {

	ComputeGS(this->_mat, this->_mu, this->_bstarabs);
	
	for (int i = 0; i < this->getNumRows(); ++i)
		this->_mu[i][i] = FT(1);
	
	this->_bstar = inv(_mu) * conv < Mat<FT> > (this->_mat);
}

/**
 *
 */

template<class T, class FT>
FT p3Matrix<T, FT>::calc_S_ik(Mat<FT>& fmat, std::vector<FT>& D_j_k, int i, int k) {
	FT ret = FT(0);

	for(int j=i; j<=k-1; j++) {
		ret += _mu[k][j] * _mu[k][j] * _bstarabs[j] * (_bstarabs[j]/D_j_k[j] - 1);
	}

	return ret;
}

/**
 *
 */
template<class T, class FT>
int p3Matrix<T, FT>::argmax_S_lk(Mat<FT>& fmat, int k, std::vector<FT>& S_lk, int& i) {
	i = -1;

	// Calculate all S_lk with 1 <= l <= k-1

	// 1. Intermediate vector with D_lk values
	std::vector<FT> D_lk;
	for(int l=0; l <= k; l++) {
		D_lk.push_back(calc_D_lk(fmat, l, k));
	}

	// 2. Fill the vector with S_lk-values
	for(int l=0; l <= k-1; l++) {
		S_lk.push_back(calc_S_ik(fmat, D_lk, l, k));
	}

	// get the index of the maximum entry
	auto max_it = std::max_element(S_lk.begin(), S_lk.end());
	i = std::distance(S_lk.begin(),max_it);

	// print for test
	for(auto it=S_lk.begin(); it!=S_lk.end(); ++it) {
		std::cout << *it << ", ";
	}
	std::cout << std::endl;

	assert (i >= 0);
	assert (S_lk.size() > 0);

	return 0;
}


/**
 * Calculates D_l^(k) as defined in [1], Page 2493, Theorem 3
 */
template<class T, class FT>
FT p3Matrix<T, FT>::calc_D_lk(Mat<FT>& fmat, int l, int k) {

	FT res = FT(0);
	for(int j=l; j<=k; j++) {
		res += (_mu[k][j] * _mu[k][j]) * (_bstarabs[j]);
	}	
	return res;
}

/**
 * Calculates SS(B) as defined in [1], Page 2490, after Theorem 1
 */
template<class T, class FT>
FT p3Matrix<T, FT>::calc_SS_B () {
	FT res = FT(1);
	for (auto it = _bstarabs.begin(); it != _bstarabs.end(); ++it) {
		res *= (*it);
	}

	return res;
}


} // end namespace p3


