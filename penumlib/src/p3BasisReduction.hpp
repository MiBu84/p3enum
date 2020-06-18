/*
 * p3BasisReduction.hpp
 *
 *  Created on: 17.06.2020
 *      Author: mburger
 */

#include "p3Matrix.hpp"

#ifndef P3BASISREDUCTION_H
#define P3BASISREDUCTION_H
namespace p3 {

template<typename T, typename FT>
struct GSInfo {
	p3Matrix<FT> _mu;
	p3Matrix<FT> _bstar;
	p3Vector<FT> _bstarabs;
	p3Vector<FT> _babs;

	GSInfo(int r, int c) {
		_mu.SetDims(r, c);
		_bstar.SetDims(r,  c);
		_bstarabs.SetLength(r);
		_babs.SetLength(r);
	}
};

template<class T, class FT>
class p3BasisReduction {
public:
	static void S2LLL (p3Matrix<T>& mat, double eta);
	static void LLL_FP (p3Matrix<T>& mat, double delta_);
	static void LLL_FP_COHEN (p3Matrix<T>& mat, FT delta_);

	static void GramSchmidt (p3Matrix<T>& mat, p3Matrix<FT>& fmat, GSInfo<T, FT>& gs_info);
	static FT calc_SS_B (const GSInfo<T, FT>& gs_info);
	static void hello(T i);

private:
	static void afterDeepInsertUpdate (const int i, const int k, const int n, GSInfo<T, FT>& gsinfo);

	static void sizeReduce(const int k, p3Matrix<T>& mat, p3Matrix<FT>& fmat, GSInfo<T, FT>& gsinfo);
	static void reduce_RET(const int k, const int l, p3Matrix<T>& mat, p3Matrix<FT>& fmat, p3Matrix<T>& H, GSInfo<T, FT>& gs_info);

	static void swapBaseVectors (p3Matrix<FT>& fmat, int idx1, int idx2);
	static void swap_SWAP(const int k, const int k_max, p3Matrix<T>& mat, p3Matrix<FT>& fmat, p3Matrix<T>& H,  GSInfo<T, FT>& gs_info);

	static int argmax_S_lk(p3Matrix<FT>& fmat, GSInfo<T, FT>& gs_info, int k, std::vector<FT>& S_lk, int& i);
	static int sigma_i_k(int i, int k, p3Matrix<T>& mat, p3Matrix<FT>& fmat);
	static FT calc_D_lk(p3Matrix<FT>& fmat, GSInfo<T,FT>& gsinfo, int l, int k);
	static FT calc_S_ik(p3Matrix<FT>& fmat, GSInfo<T, FT>& gs_info, std::vector<FT>& D_j_k, int i, int k);
};

#include "p3BasisReduction.tcc"

} // namespace p3
#endif





