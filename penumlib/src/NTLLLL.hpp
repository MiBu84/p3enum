
#ifndef NTLLLL_H
#define NTLLLL_H

#include <NTL/mat_ZZ.h>

//NTL_START_IMPL

using namespace NTL;

class NTLLLL {
public:
	static long LLL(ZZ& det, mat_ZZ& B, long a, long b, long verbose);
	static long LLL (vec_ZZ& D, mat_ZZ& B, mat_ZZ* U, long a, long b, long verbose);
	static void IncrementalGS(mat_ZZ& B, vec_long& P, vec_ZZ& D, vec_vec_ZZ& lam, long& s, long k);

private:
	static void reduce(long k, long l,
	            mat_ZZ& B, vec_long& P, vec_ZZ& D,
	            vec_vec_ZZ& lam, mat_ZZ* U);

	static long SwapTest(const ZZ& d0, const ZZ& d1, const ZZ& d2, const ZZ& lam,
	                     long a, long b);

	static long swap(long k, mat_ZZ& B, vec_long& P, vec_ZZ& D,
	          vec_vec_ZZ& lam, mat_ZZ* U, long m, long verbose);

	static void RowTransform(ZZ& c1, ZZ& c2,
	                         const ZZ& x, const ZZ& y, const ZZ& u, const ZZ& v);
	static void RowTransform(vec_ZZ& c1, vec_ZZ& c2,
	                         const ZZ& x, const ZZ& y, const ZZ& u, const ZZ& v);

	static void BalDiv(ZZ& q, const ZZ& a, const ZZ& d);
	static void MulSubFrom(vec_ZZ& c, const vec_ZZ& c2, const ZZ& x);
	static void MulSubFrom(vec_ZZ& c, const vec_ZZ& c2, long x);
	static void MulAddDiv(ZZ& c, const ZZ& c1, const ZZ& c2, const ZZ& x, const ZZ& y, const ZZ& z);
	static void MulSubDiv(ZZ& c, const ZZ& c1, const ZZ& c2, const ZZ& x, const ZZ& y, const ZZ& z);
	static void ExactDiv(ZZ& qq, const ZZ& a, const ZZ& b);

};

#endif

//NTL_END_IMPL
