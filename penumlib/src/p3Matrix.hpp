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
#include "p3Vector.hpp"

using namespace NTL;

namespace p3 {
template<class T>
	class p3Matrix {
	public:
		p3Matrix();
		p3Matrix(int r, int c, bool is_ident=false);
		p3Matrix(const NTL::Mat<ZZ>& mat);

		p3Matrix& operator=(const p3Matrix& other) {
			this->_data = other._data;
			return *this;
		}

		p3Vector<T>& operator[] (int idx) {
			return _data[idx];
		}

		friend std::ostream& operator<<(ostream& os, p3Matrix<T>& dt) {
			os << "[ ";

			for (unsigned int i=0; i<dt.NumRows(); i++) {
				os << dt[i] << " ";
			}

			os  << endl << "]" << endl;

			return os;
		}

		unsigned int getNumRows() {
			return _data.size();
		}

		unsigned int getNumCols() {
			if (_data.size() > 0)
				return _data[0].size();
			return 0;
		}

		unsigned int NumRows() {
			return _data.size();
		}

		unsigned int NumCols() {
			if (_data.size() > 0)
				return _data[0].size();
			return 0;
		}

		int getDim() {
			return getNumCols();
		}

		int getRank() {
			return getNumRows();
		}

		void resize(int r, int c) {
			_data.resize(r);
			for(auto it=_data.begin(); it!=_data.end(); ++it) {
				it->resize(c);
			}
		}

		void SetDims(int r, int c) {
			resize(r,c);
		}

		//void BKZ();

		//void GramSchmidt();
		//void ModifiedGramSchmidt();
		//void GramSchmidtIntegral();
		//void afterDeepInsertUpdate(int i, int k);
		//void S2LLL(double eta);
		//void LLL_FP(double eta);

		//void GSO();
		//void GSO_FP();

		//p3Matrix<T> getNTLMat() {
		//	return _mat;
		//}


	private:
		/*void sizeReduce(int const k, p3Matrix<FT>& fmat);

		int sigma_i_k(int i, int k, p3Matrix<FT>& fmat);
		int argmax_S_lk(p3Matrix<FT>& fmat, int k, std::vector<FT>& S_i_k, int& i);

		FT calc_D_lk(p3Matrix<FT>& fmat, int l, int k);
		FT calc_S_ik(p3Matrix<FT>& fmat, std::vector<FT>& S_i_k, int i, int k);
		FT calc_SS_B();

		void swapBaseVectors (int idx1, int idx2);*/

		std::vector<p3Vector<T>> _data;

		//p3Matrix<T> _mat;;
		//p3Matrix<FT> _fmat; // The floating point representation of the matrix

		//p3Matrix<FT> _mu;
		//p3Matrix<FT> _bstar;
		//Vec<FT> _bstarabs;
	};
}

#include "p3Matrix.tcc"

#endif /* SRC_P3MATRIX_HPP_ */
