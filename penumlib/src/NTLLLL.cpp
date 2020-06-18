/*
 * NTLLLL.cpp
 *
 *  Created on: 16.06.2020
 *      Author: mburger
 */

#include "NTLLLL.hpp"

#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
#include <iostream>

using namespace NTL;
using namespace std;

long NTLLLL::SwapTest(const ZZ& d0, const ZZ& d1, const ZZ& d2, const ZZ& lam,  long a, long b)

// test if a*d1^2 > b*(d0*d2 + lam^2)

{
   NTL_ZZRegister(t1);
   NTL_ZZRegister(t2);

   mul(t1, d0, d2);
   sqr(t2, lam);
   add(t1, t1, t2);
   mul(t1, t1, b);

   sqr(t2, d1);
   mul(t2, t2, a);

   return t2 > t1;
}

void NTLLLL::ExactDiv(ZZ& qq, const ZZ& a, const ZZ& b)
{
   NTL_ZZRegister(q);
   NTL_ZZRegister(r);

   DivRem(q, r, a, b);
   if (!IsZero(r)) {
      cerr << "a = " << a << "\n";
      cerr << "b = " << b << "\n";
      LogicError("ExactDiv: nonzero remainder");
   }
   qq = q;
}

void NTLLLL::MulSubDiv(ZZ& c, const ZZ& c1, const ZZ& c2,
                      const ZZ& x, const ZZ& y, const ZZ& z)

// c = (x*c1 - y*c2)/z

{
   NTL_ZZRegister(t1);
   NTL_ZZRegister(t2);

   mul(t1, x, c1);
   mul(t2, y, c2);
   sub(t1, t1, t2);
   ExactDiv(c, t1, z);
}

void NTLLLL::MulAddDiv(ZZ& c, const ZZ& c1, const ZZ& c2,
                      const ZZ& x, const ZZ& y, const ZZ& z)

// c = (x*c1 + y*c2)/z

{
   NTL_ZZRegister(t1);
   NTL_ZZRegister(t2);

   mul(t1, x, c1);
   mul(t2, y, c2);
   add(t1, t1, t2);
   ExactDiv(c, t1, z);
}

void NTLLLL::MulSubFrom(vec_ZZ& c, const vec_ZZ& c2, long x)

// c = c - x*c2

{
   long n = c.length();
   if (c2.length() != n) LogicError("MulSubFrom: length mismatch");

   long i;
   for (i = 1; i <= n; i++)
      NTL::MulSubFrom(c(i), c2(i), x);
}

void NTLLLL::MulSubFrom(vec_ZZ& c, const vec_ZZ& c2, const ZZ& x)
// c = c - x*c2

{
   long n = c.length();
   if (c2.length() != n) LogicError("MulSubFrom: length mismatch");

   long i;
   for (i = 1; i <= n; i++)
	   NTL::MulSubFrom(c(i), c2(i), x);
}

void NTLLLL::BalDiv(ZZ& q, const ZZ& a, const ZZ& d)

//  rounds a/d to nearest integer, breaking ties
//    by rounding towards zero.  Assumes d > 0.

{
   NTL_ZZRegister(r);
   DivRem(q, r, a, d);


   add(r, r, r);

   long cmp = compare(r, d);
   if (cmp > 0 || (cmp == 0 && q < 0))
      add(q, q, 1);
}

void NTLLLL::RowTransform(vec_ZZ& c1, vec_ZZ& c2,
                         const ZZ& x, const ZZ& y, const ZZ& u, const ZZ& v)

// (c1, c2) = (x*c1 + y*c2, u*c1 + v*c2)

{
   long n = c1.length();
   if (c2.length() != n) LogicError("MulSubDiv: length mismatch");
   NTL_ZZRegister(t1);
   NTL_ZZRegister(t2);
   NTL_ZZRegister(t3);
   NTL_ZZRegister(t4);

   long i;
   for (i = 1; i <= n; i++) {
      mul(t1, x, c1(i));
      mul(t2, y, c2(i));
      add(t1, t1, t2);

      mul(t3, u, c1(i));
      mul(t4, v, c2(i));
      add(t3, t3, t4);

      c1(i) = t1;
      c2(i) = t3;
   }
}

void NTLLLL::RowTransform(ZZ& c1, ZZ& c2,
                         const ZZ& x, const ZZ& y, const ZZ& u, const ZZ& v)

// (c1, c2) = (x*c1 + y*c2, u*c1 + v*c2)

{
   NTL_ZZRegister(t1);
   NTL_ZZRegister(t2);
   NTL_ZZRegister(t3);
   NTL_ZZRegister(t4);

   mul(t1, x, c1);
   mul(t2, y, c2);
   add(t1, t1, t2);

   mul(t3, u, c1);
   mul(t4, v, c2);
   add(t3, t3, t4);

   c1 = t1;
   c2 = t3;
}

long NTLLLL::swap(long k, mat_ZZ& B, vec_long& P, vec_ZZ& D,
          vec_vec_ZZ& lam, mat_ZZ* U, long m, long verbose)

// swaps vectors k-1 and k;  assumes P(k-1) != 0
// returns 1 if vector k-1 need to be reduced after the swap...
//    this only occurs in 'case 2' when there are linear dependencies

{
   long i, j;
   NTL_ZZRegister(t1);
   NTL_ZZRegister(t2);
   NTL_ZZRegister(t3);
   NTL_ZZRegister(e);
   NTL_ZZRegister(x);
   NTL_ZZRegister(y);


   if (P(k) != 0) {
      if (verbose) cerr << "swap case 1: " << k << "\n";

      NTL::swap(B(k-1), B(k));
      if (U) NTL::swap((*U)(k-1), (*U)(k));

      for (j = 1; j <= k-2; j++)
         if (P(j) != 0)
            NTL::swap(lam(k-1)(P(j)), lam(k)(P(j)));

      for (i = k+1; i <= m; i++) {
         MulAddDiv(t1, lam(i)(P(k)-1), lam(i)(P(k)),
                   lam(k)(P(k)-1), D[P(k)-2], D[P(k)-1]);
         MulSubDiv(t2, lam(i)(P(k)-1), lam(i)(P(k)),
                   D[P(k)], lam(k)(P(k)-1), D[P(k)-1]);
         lam(i)(P(k)-1) = t1;
         lam(i)(P(k)) = t2;
      }

      MulAddDiv(D[P(k)-1], D[P(k)], lam(k)(P(k)-1),
                D[P(k)-2], lam(k)(P(k)-1), D[P(k)-1]);

      return 0;
   }
   else if (!IsZero(lam(k)(P(k-1)))) {
      if (verbose) cerr << "swap case 2: " << k << "\n";
      XGCD(e, x, y, lam(k)(P(k-1)), D[P(k-1)]);

      ExactDiv(t1, lam(k)(P(k-1)), e);
      ExactDiv(t2, D[P(k-1)], e);

      t3 = t2;
      NTL::negate(t2, t2);
      RowTransform(B(k-1), B(k), t1, t2, y, x);
      if (U) RowTransform((*U)(k-1), (*U)(k), t1, t2, y, x);
      for (j = 1; j <= k-2; j++)
         if (P(j) != 0)
            RowTransform(lam(k-1)(P(j)), lam(k)(P(j)), t1, t2, y, x);

      sqr(t2, t2);
      ExactDiv(D[P(k-1)], D[P(k-1)], t2);

      for (i = k+1; i <= m; i++)
         if (P(i) != 0) {
            ExactDiv(D[P(i)], D[P(i)], t2);
            for (j = i+1; j <= m; j++) {
               ExactDiv(lam(j)(P(i)), lam(j)(P(i)), t2);
            }
         }

      for (i = k+1; i <= m; i++) {
         ExactDiv(lam(i)(P(k-1)), lam(i)(P(k-1)), t3);
      }

      NTL::swap(P(k-1), P(k));

      return 1;
   }
   else {
      if (verbose) cerr << "swap case 3: " << k << "\n";

      NTL::swap(B(k-1), B(k));
      if (U) NTL::swap((*U)(k-1), (*U)(k));

      for (j = 1; j <= k-2; j++)
         if (P(j) != 0)
            NTL::swap(lam(k-1)(P(j)), lam(k)(P(j)));

      NTL::swap(P(k-1), P(k));

      return 0;
   }
}


void NTLLLL::reduce(long k, long l,
            mat_ZZ& B, vec_long& P, vec_ZZ& D,
            vec_vec_ZZ& lam, mat_ZZ* U)
{
   NTL_ZZRegister(t1);
   NTL_ZZRegister(r);

   if (P(l) == 0) return;
   add(t1, lam(k)(P(l)), lam(k)(P(l)));
   abs(t1, t1);
   if (t1 <= D[P(l)]) return;

   long j;
   long rr, small_r;

   BalDiv(r, lam(k)(P(l)), D[P(l)]);

   if (r.WideSinglePrecision()) {
      small_r = 1;
      rr = to_long(r);
   }
   else {
      small_r = 0;
   }

   if (small_r) {
      MulSubFrom(B(k), B(l), rr);

      if (U) MulSubFrom((*U)(k), (*U)(l), rr);

      for (j = 1; j <= l-1; j++)
         if (P(j) != 0)
            NTL::MulSubFrom(lam(k)(P(j)), lam(l)(P(j)), rr);
      NTL::MulSubFrom(lam(k)(P(l)), D[P(l)], rr);
   }
   else {
      MulSubFrom(B(k), B(l), r);

      if (U) MulSubFrom((*U)(k), (*U)(l), r);

      for (j = 1; j <= l-1; j++)
         if (P(j) != 0)
        	 NTL::MulSubFrom(lam(k)(P(j)), lam(l)(P(j)), r);
      NTL::MulSubFrom(lam(k)(P(l)), D[P(l)], r);
   }


}

void NTLLLL::IncrementalGS(mat_ZZ& B, vec_long& P, vec_ZZ& D, vec_vec_ZZ& lam, long& s, long k) {
	   //long n = B.NumCols();
	   //long m = B.NumRows();

	   NTL_ZZRegister(u);
	   NTL_ZZRegister(t1);
	   NTL_ZZRegister(t2);

	   long i, j;

	   for (j = 1; j <= k-1; j++) {
	      long posj = P(j);
	      if (posj == 0) continue;

	      InnerProduct(u, B(k), B(j));
	      for (i = 1; i <= posj-1; i++) {
	         mul(t1, D[i], u);
	         mul(t2, lam(k)(i), lam(j)(i));
	         sub(t1, t1, t2);
	         div(t1, t1, D[i-1]);
	         u = t1;
	      }

	      lam(k)(posj) = u;
	   }

	   InnerProduct(u, B(k), B(k));
	   for (i = 1; i <= s; i++) {
	      mul(t1, D[i], u);
	      mul(t2, lam(k)(i), lam(k)(i));
	      sub(t1, t1, t2);
	      div(t1, t1, D[i-1]);
	      u = t1;
	   }

	   if (u == 0) {
	      P(k) = 0;
	   }
	   else {
	      s++;
	      P(k) = s;
	      D[s] = u;
	   }
}


long NTLLLL::LLL (vec_ZZ& D, mat_ZZ& B, mat_ZZ* U, long a, long b, long verbose) {
   long m = B.NumRows();
   //long n = B.NumCols();

   long force_reduce = 1;

   vec_long P;
   P.SetLength(m);

   D.SetLength(m+1);
   D[0] = 1;

   vec_vec_ZZ lam;

   lam.SetLength(m);

   long j;
   for (j = 1; j <= m; j++)
      lam(j).SetLength(m);

   if (U) ident(*U, m);

   long s = 0;

   long k = 1;
   long max_k = 0;


   while (k <= m) {
      if (k > max_k) {
         IncrementalGS(B, P, D, lam, s, k);
         max_k = k;
      }

      if (k == 1) {
         force_reduce = 1;
         k++;
         continue;
      }

      if (force_reduce)
         for (j = k-1; j >= 1; j--)
            reduce(k, j, B, P, D, lam, U);

      if (P(k-1) != 0 &&
          (P(k) == 0 ||
           SwapTest(D[P(k)], D[P(k)-1], D[P(k)-2], lam(k)(P(k)-1), a, b))) {
         force_reduce = swap(k, B, P, D, lam, U, max_k, verbose);
         k--;
      }
      else {
         force_reduce = 1;
         k++;
      }
      //cout << "Wichsdoedel" << endl;
      //break;
   }

   D.SetLength(s+1);
   return s;
}

long NTLLLL::LLL(ZZ& det, mat_ZZ& B, long a, long b, long verbose)
{
   if (a <= 0 || b <= 0 || a > b || b/4 >= a) LogicError("LLL: bad args");

   vec_ZZ D;
   long s;
   s = LLL(D, B, 0, a, b, verbose);
   det = D[s];
   return s;
}

