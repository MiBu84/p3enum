/*
 * p3Matrix.hpp
 *
 *  Created on: 08.06.2020
 *      Author: mburger
 */

#ifndef SRC_P3MATRIX_HPP_
#define SRC_P3MATRIX_HPP_

#include <NTL/mat_ZZ.h>
#include <NTL/matrix.h>
#include <NTL/vector.h>
#include <vector>

using namespace NTL;

namespace p3 {
template<class T, class FT>
	class p3Matrix {
	public:
		p3Matrix();
		p3Matrix(int r, int c);
		p3Matrix(const Mat<T>& mat);

		p3Matrix& operator=(const p3Matrix& other) {
			this->_mat = other._mat;
			return *this;
		}

		Vec<T>& operator[] (int idx) {
			return _mat[idx];
		}

		int getNumRows() {
			return _mat.NumRows();
		}

		int getNumCols() {
			return _mat.NumCols();
		}

		int getDim() {
			return _mat.NumCols();
		}

		int getRank() {
			return _mat.NumRows();
		}

		void BKZ();

		void S2LLL(double eta);
		void LLL_FP(double eta);

		void GSO();
		void GSO_FP();


	private:
		void sizeReduce(int k, Mat<FT>& fmat);

		int sigma_i_k(int i, int k, Mat<FT>& fmat);
		int argmax_S_lk(Mat<FT>& fmat, int k, std::vector<FT>& S_i_k, int& i);

		FT calc_D_lk(Mat<FT>& fmat, int l, int k);
		FT calc_S_ik(Mat<FT>& fmat, std::vector<FT>& S_i_k, int i, int k);
		FT calc_SS_B();

		void swapBaseVectors (int idx1, int idx2);

		Mat<T> _mat;
		Mat<FT> _fmat; // The floating point representation of the matrix

		Mat<FT> _mu;
		Mat<FT> _bstar;
		Vec<FT> _bstarabs;
	};
}

#include "p3Matrix.tcc"

#endif /* SRC_P3MATRIX_HPP_ */
