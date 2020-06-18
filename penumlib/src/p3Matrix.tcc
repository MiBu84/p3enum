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

/*template<class T>
T operator* (const Vec<T>& a, const Vec<T>& b) {
	assert(a.length() == b.length());

	T res = T(0);

	auto itv1 = a.begin();
	auto itv2 = b.begin();

	for(;itv1 != a.end(); ++itv1, ++itv2) {
		res += (*itv1) * (*itv2);
	}

	return res;

}

template<class T>
Vec<T> operator* (const T a, const Vec<T>& b) {
	Vec<T> res; res.SetLength(b.length());

	auto it_res = res.begin();

	for(auto itv1 = b.begin(); itv1 != b.end(); ++itv1, ++it_res) {
		*it_res = a * (*itv1);
	}

	return res;
}

template<class T>
Vec<T> operator- (const Vec<T>& a, const Vec<T>& b) {
	assert(a.length() == b.length());
	Vec<T> res; res.SetLength(b.length());
	auto it_res = res.begin();

	auto itv1 = a.begin();
	auto itv2 = b.begin();

	for(;itv1 != a.end(); ++itv1, ++itv2, ++it_res) {
		*it_res = (*itv1) - (*itv2);
	}

	return res;

}*/

/*template <class T>
T round (const T& val) {
	if constexpr (std::is_integral_v<T>) {  // constexpr only necessary on first statement
		return T(val);
	}
	else if (std::is_floating_point_v<T>) {  // automatically constexpr {
		return std::round(val);
	}
	else {
		return T(val);
	}
}*/

/*long round (const long& val) {
	return long(val);
}*/

template<class T1, class T2>
	T2 conv_to_ft (const T1& input) {
	T2 ret = T2(0);

	std::string strval = convert2Str <T1> ( input );
	convertFromStr<T2>(ret, strval);
	return ret;
}

template<class T1, class T2>
	T2 conv_to_int (const T1& input) {
	T2 ret = T2(0);

	std::string strval = convert2Str <T1> ( (input) );
	convertFromStr<T2>(ret, strval);
	return ret;
}

template<class T1, class T2>
	T2 conv_to_int_from_ft (const T1& input) {
	T2 ret = T2(0);

	std::string strval = convert2Str <T1> ( round(input) );
	convertFromStr<T2>(ret, strval);
	return ret;
}

template<class T1, class T2>
	p3Matrix<T2> conv_to_int_mat (p3Matrix<T1>& input) {
	p3Matrix<T2> retmat;
	retmat.SetDims(input.NumRows(), input.NumCols());

	// Copy results to the integer-representation
	assert(retmat.NumRows() == input.NumRows());
	assert(retmat.NumCols() == input.NumCols());

	for(unsigned int r = 0; r<input.NumRows();r++) {
		for(unsigned int c = 0; c < input.NumCols(); c++) {
			retmat[r][c] = conv_to_int_from_ft<T1, T2>(input[r][c]);
		}
	}

	return retmat;
}

template<class T1, class T2>
	p3Matrix<T2> conv_to_ft_mat (p3Matrix<T1>& input) {

	p3Matrix<T2> retmat;
	retmat.SetDims(input.NumRows(), input.NumCols());

	// Copy results to the integer-representation
	assert(retmat.NumRows() == input.NumRows());
	assert(retmat.NumCols() == input.NumCols());

	for(unsigned int r = 0; r<input.NumRows();r++) {
		for(unsigned int c = 0; c < input.NumCols(); c++) {
			retmat[r][c] = conv_to_ft<T1, T2>(input[r][c]);
		}
	}

	return retmat;
}

template<class T>
	p3Matrix<T>::p3Matrix() {

	}

template<class T>
	p3Matrix<T>::p3Matrix(int r, int c, bool is_ident) {
		resize(r, c);

		for(unsigned int r=0; r < this->NumRows(); r++) {
			for(unsigned  int c=0; c < this->NumRows(); c++) {
				this->_data[r][c] = T(0);
				// If user wants an identity matrix
				if (is_ident && r==c) {
					this->_data[r][c] = T(1);
				}
			}
		}
	}

template<class T>
	p3Matrix<T>::p3Matrix (const NTL::Mat<ZZ>& mat) {

		resize(mat.NumRows(), mat.NumCols());

		for(int r=0; r < mat.NumRows(); r++) {
			for(int c=0; c < mat.NumRows(); c++) {
				this->_data[r][c] = p3::conv_to_int<ZZ, T> (mat[r][c]);
			}
		}

		//this->_bstar.SetDims(mat.NumRows(), mat.NumCols());
		//this->_mu.SetDims(mat.NumRows(), mat.NumCols());
		//this->_bstarabs.SetLength(mat.NumCols());
	}

/*template<class T, class FT>
	void p3Matrix<T, FT>::BKZ () {	
}*/

/*template<class T, class FT>
	void p3Matrix<T, FT>::ModifiedGramSchmidt () {

	int n = _mat.NumRows();
	int m = _mat.NumCols();

	std::vector<mat_ZZ> Aps;

	mat_ZZ lmu; lmu.SetDims(n, m);

	for(int k=0;  k < n; k++) {
		Aps.push_back(_mat);
	}

	for (int k = 0; k < n; k++) {
		lmu[k][k] = Aps[k][k] * Aps[k][k];
	}

}*/

/*template<class T, class FT>
	void p3Matrix<T, FT>::afterDeepInsertUpdate (const int i, const int k) {

	int n = _fmat.NumRows();
	FT tmp1;

	Vec<FT> D, Dinv, P, epsilons, S;
	D.SetLength(n); Dinv.SetLength(n);
	P.SetLength(n);
	epsilons.SetLength(n);
	S.SetLength(n);

	P[k] = _bstarabs[k]; // P_k = B_k
	D[k] = _bstarabs[k]; // D_k = B_k

	// 2:4
	for (int j = k - 1; j >= i; --j) {
		P[j] = _mu[k][j] * _bstarabs[j];
		D[j] = D[j + 1] + _mu[k][j] * P[j];
	}


	// 5 Set S_i = S_i+1 = ... = S_n = 0
	for (int c = i; c < n; c++) {
		S[c] = 0;
	}

	// 30:32
	Dinv[k] = 1.0 / D[k]; // da mu[k][k] == 1

	for (int j = k; j >= i + 1; --j) {
		Dinv[j-1] = 1.0 / D[j-1];
		_bstarabs[j] = (D[j] * _bstarabs[j-1]) * Dinv[j-1];
	}
	_bstarabs[i] = D[i];

	// 8
	tmp1 = _mu[k][k-1] * Dinv[k];
	// 9:11
	for (int l = n - 1; l >= k+1; --l) {
		S[l] = P[k] * _mu[l][k];
		_mu[l][k] = _mu[l][k-1] - tmp1 * S[l];
	}

	// 7
	for (int j = k - 1; j >= i+1; --j) {
		// 8
		tmp1 = _mu[k][j-1] * Dinv[j];

		// 9:11
		for (int l = k + 1; l < n; ++l) {
			S[l] += P[j] * _mu[l][j];
			_mu[l][j] = _mu[l][j-1] - tmp1 * S(l);
		}

		// 12:14
		for (int l = k; l >= j + 2; --l) {
			S(l) += P(l) * _mu[l-1][j];
			_mu[l][j] = _mu[l-1][j-1] - tmp1 * S(l);
		}

		// 15:17
		S(j+1) = P(j);
		_mu[j+1][j] = _mu[j][j-1] - tmp1 * S[j+1];
	}

	// 21
	for (int l = k + 1; l < n; ++l) {
		_mu[l][i] = (S[l] + P[i] * _mu[l][i]) * Dinv[i];
	}

	// 22
	for (int l = k; l >= i+2; --l) {
		_mu[l][i] = (S[l] + P[i] * _mu[l-1][i]) * Dinv[i];
	}

	// 23
	_mu[i+1][i] = P[i] * Dinv[i];

	// Copy von Zeile 26
	for(int j = 1; j <= i-1; ++j) {
		epsilons[j] = _mu[k][j];
	}

	// 25:29
	for(int j = k; j != i; --j) {
		for(int l = 1; l < i; ++l) {
			_mu[j][l] = _mu[j-1][l];
		}
	}

	// Copy von Zeile 28
	for(int j = 1; j < i; ++j) {
		_mu[i][j] = epsilons[j];
	}

}*/

/*template<class T, class FT>
	void p3Matrix<T, FT>::GramSchmidt () {

	this->_bstar[0] = this->_fmat[0];
	for(int i=1; i < this->_fmat.NumRows(); i++) {
		Vec<FT> v = _fmat[i];

		for(int j=i-1; j>= 0; j--) {
			_mu[i][j] = (_fmat[i] * _bstar[j]) / ( _bstar[j] * _bstar[j]);
			v = v - (_mu[i][j] * _bstar[j]);
		}
		_bstar[i] = v;
	}

	for(int i=0; i<this->_fmat.NumRows(); i++) {
		_mu[i][i] = 1.0;
	}

	// Fill ||b*||^2
	assert(_bstarabs.length() == _bstar.NumCols());
	auto it_abs=_bstarabs.begin();
	for(int c=0; c < _bstar.NumCols(); c++, ++it_abs) {
		*it_abs = _bstar[c] * _bstar[c];
	}
	//cout << _mu << endl;

	return;
}*/

/*template<class T, class FT>
	void p3Matrix<T, FT>::LLL_FP (double delta_) {
	FT delta = FT(delta_);

	_fmat = conv_to_ft_mat <T, FT> (_mat);
	GramSchmidt();
	int k=1;
	int redcnt = 0;
	FT startBB = calc_SS_B();
	cout << "Start SS " << startBB << endl;
	while (k < getDim()) {
		sizeReduce(k, _fmat);
		redcnt++;

		if (_bstarabs[k] >= (delta - _mu[k][k-1] * _mu[k][k-1]) * _bstarabs[k-1] ) {
			k = k + 1;
		}

		else {
			swapBaseVectors(k, k-1);
			GramSchmidt();
			k = std::max(1,k-1);
		}
	}

	cout << _fmat << endl;
	FT endBB = calc_SS_B();
	std::cout << "Improved to: " << endBB << " / " << endBB/startBB << std::endl;
	std::cout << redcnt << " reductions." << std::endl;
}*/

/*template<class T, class FT>
	void p3Matrix<T, FT>::S2LLL (double eta) {
	_fmat = conv_to_ft_mat <T, FT> (_mat);
	GramSchmidt();

	FT startBB = calc_SS_B();
	cout << startBB << endl;
	int redcnt = 0;

	int k = 1;
	while (k < getDim() ) {

		//std::cout << "Old: " << _mu << endl;
		sizeReduce(k, _fmat);
		//std::cout << "New: " << _mu << endl;
		redcnt++;

		// Compute i <-argmax
		int i = -1;
		std::vector<FT> S_jk;//NTL::conv<FT>(std::numeric_limits<double>min());
		argmax_S_lk(_fmat, k, S_jk, i);

		if (S_jk[i] <= (1-eta) * calc_SS_B()) {
			k++;
		}

		else {
			//cout << "Deep with " << i << " / " << k << endl;
			sigma_i_k(i,k,_fmat);
			afterDeepInsertUpdate(i, k);
			k = std::max(i, 1);
		}
	}

	// Copy results to the integer-representation
	_mat = conv_to_int_p3Matrix<FT, T> (_fmat);

	std::cout << _fmat << endl;
	FT endBB = calc_SS_B();
	std::cout << "Improved to: " << endBB << " / " << endBB/startBB << std::endl;
	std::cout << redcnt << " reductions." << std::endl;
}*/

/*template<class T, class FT>
	void p3Matrix<T, FT>::sizeReduce(const int k, p3Matrix<FT>& fmat) {

	//cout << "k: " << k << " / ";
	for(int j =  k-1; j >= 0; j--) {

		if(fabs(_mu[k][j]) > 0.5) {
			FT q_j = round( _mu[k][j] );
			//cout << j << ",";

			// RR to ZZ
			//std::string strval = convert2Str<FT>(q_j);
			//T q_j_i; convertFromStr<T>(q_j_i, strval);
			//_mat[k] = _mat[k] - q_j_i * _mat[j];

			fmat[k] = fmat[k] - q_j * fmat[j];
			//_fmat = conv_to_ft_mat <T, FT> (_mat);

			// Update only relevant mu
		    for (int i = 0; i < j; i++)
		      _mu[k][i] = _mu[k][i] - q_j * _mu[j][i];
			//GramSchmidt();
		}
		// Full update
		GramSchmidt();
	}
	//cout << endl;
}*/

/*template<class T, class FT>
void p3Matrix<T, FT>::swapBaseVectors (int idx1, int idx2) {
	Vec<FT> vec1 = _fmat[idx1];
	Vec<FT> vec2 = _fmat[idx2];
	_fmat[idx1] = vec2;
	_fmat[idx2] = vec1;

	//Vec<T> vec1i = _mat[idx1];
	//Vec<T> vec2i = _mat[idx2];
	//_mat[idx1] = vec2i;
	//_mat[idx2] = vec1i;

	return;
}*/

/*template<class T, class FT>
int p3Matrix<T, FT>::sigma_i_k(int i, int k, p3Matrix<FT>& fmat) {

	//if(abs(i - k ) > 1) {
		//std::cout << "Doing it deep " << k-i << std::endl;
	//}

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
}*/

/*template<class T, class FT>
	void p3Matrix<T, FT>::GSO() {

	ComputeGS(this->_mat, this->_mu, this->_bstarabs);
	

	for (int i = 0; i < this->getNumRows(); ++i)
		this->_mu[i][i] = FT(1);
	
	this->_bstar = inv(_mu) * conv < p3Matrix<FT> > (this->_mat);

	//cout << this->_bstarabs << endl;
}*/

/**
 *
 */

/*template<class T, class FT>
FT p3Matrix<T, FT>::calc_S_ik(p3Matrix<FT>& fmat, std::vector<FT>& D_j_k, int i, int k) {
	FT ret = FT(0);

	for(int j=i; j<=k-1; j++) {
		ret += _mu[k][j] * _mu[k][j] * _bstarabs[j] * (_bstarabs[j]/D_j_k[j] - 1);
	}

	return ret;
}*/

/**
 *
 */
/*template<class T, class FT>
int p3Matrix<T, FT>::argmax_S_lk(p3Matrix<FT>& fmat, int k, std::vector<FT>& S_lk, int& i) {
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

	assert (i >= 0);
	assert (S_lk.size() > 0);

	return 0;
}*/


/**
 * Calculates D_l^(k) as defined in [1], Page 2493, Theorem 3
 */
/*template<class T, class FT>
FT p3Matrix<T, FT>::calc_D_lk(p3Matrix<FT>& fmat, int l, int k) {

	FT res = FT(0);
	for(int j=l; j<=k; j++) {
		res += (_mu[k][j] * _mu[k][j]) * (_bstarabs[j]);
	}	
	return res;
}*/

/**
 * Calculates SS(B) as defined in [1], Page 2490, after Theorem 1
 */
/*template<class T, class FT>
FT p3Matrix<T, FT>::calc_SS_B () {
	FT res = FT(0);
	for (auto it = _bstarabs.begin(); it != _bstarabs.end(); ++it) {
		res += (*it);
	}

	return res;
}*/


} // end namespace p3


