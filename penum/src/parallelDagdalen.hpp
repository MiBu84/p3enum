/*
 * parallelDagdalen.hpp
 *
 *  Created on: 15.05.2018
 *      Author: tuxed
 */

#ifndef SRC_PARALLELDAGDALEN_HPP_
#define SRC_PARALLELDAGDALEN_HPP_

#include <NTL/RR.h>
#include <NTL/mat_RR.h>
#include <NTL/vector.h>
#include <vector>

using namespace std;
using namespace NTL;

template <class ZT, class FT>
 struct subEnumDef {
	Vec<ZZ> ubar;
	FT lbar;
	FT cbar;
	int tbar;
};

RR glob_A;
Vec<ZZ> glob_umin;

vector<subEnumDef<ZZ, RR>> L;







void printBasis(const Mat<ZZ>& B) {
	cout << "[ " << endl;
	for(int r=0; r<B.NumRows(); r++) {
		for(int c=0; c<B.NumCols(); c++) {
			cout << B[r][c] << "\t";
		}
		cout << endl;
	}
	cout << "]" << endl;
}

template<class ZT, class FT> void DagdelenSubEnumeration(const Mat<ZT>& B, const Mat<FT>& mu,
			const Vec<FT>& bstarnorm, const FT& A_in, int j, int k,
			int dim, const subEnumDef<ZT, FT>& sub_info) {

	// Set scalars
	int bound = k;
	int s, t;
	t = s = sub_info.tbar;

	// Vector types
	Vec<ZT> u;	u.SetLength(k + 2);
	u = sub_info.ubar;

	//Vec<ZT> umin;	umin.SetLength(k + 2);
	//umin = u;

	Vec<FT> l; l.SetLength(k + 2);
	Vec<FT> c; c.SetLength(k + 2);

	// To decide the zigzag pattern
	Vec<ZT> Delta; Delta.SetLength(k + 2);
	Vec<ZT> delta; delta.SetLength(k + 2);
	Vec<ZT> v; v.SetLength(k + 2);

	for(int i=0; i <= k; i++) {
		c[i] = l[i] = 0;
		Delta[i] = v[t] = 0;
		delta[i] = 1;
	}

	// Set vectors
	l[t + 1] = sub_info.lbar;
	c[t] = sub_info.cbar;

	v[j] = 0;
	Delta[j] = 0;
	delta[j] = 1;

	for(int i=j+1; i <= k; i++) {
		c[i] = l[i] = 0;
		u[i] = Delta[i] = v[t] = 0;
		delta[i] = 1;
	}

	//s = t = j;
	//RR A = 1.05 * bstarnorm[0];

	while(t <= bound) {
		l[t] = l[ t + 1 ] + ( conv<FT>(u[t]) + c[t] ) * ( conv<FT>(u[t]) + c[t] ) * bstarnorm[t];
		cout <<"L: " << l << " / U: " << u << " / C: " << c << " / t: "<< t << endl;
		if( l[t] < glob_A) {
			if (t > j) {
				t = t - 1;
				c[t] = muProd<ZT, FT>(u, mu, t, k);
				u[t] = v[t] = conv<ZT>(ceil(-c[t] - 0.5));

				// For zigzag
				Delta[t] = 0;
				if(conv<FT>(u[t]) > -c[t]) {
					delta[t] = - 1;
				}
				else {
					delta[t] = 1;
				}

				if(bound == k) {
					subEnumDef<ZT, FT> new_def;
					new_def.ubar = u;
					new_def.lbar = l[t+2];
					new_def.cbar = c[t+1];
					new_def.tbar = t + 1;
					L.push_back(new_def);
					cout << "Pushing treeInfo"
							<< "((" << u << "), "
							<< t + 1 << ")" << endl;
					bound = t;
				}
			}

			else {
				glob_A = l[t];
				glob_umin = u;
			}
		}

		else {
			t = t + 1;

			// choose next value for u t using the zig-zag pattern
			s = max(s,t);
			if(t<s)
				Delta[t] = -Delta[t];
			if(Delta[t] * delta[t] >= 0)
				Delta[t] += delta[t];
			u[t] = v[t] + Delta[t];
		}
	}
	//calcVectorLength<ZT,FT>(umin, B);
	return;
}


template<class ZT, class FT> void DagdelenEnumerationMaster(const Mat<ZZ>& B, const Mat<RR>& mu,
			const Vec<FT>& bstarnorm, const FT& A_in, int j, int k, int dim) {

	Vec<ZT> u; u.SetLength(k - j + 2);
	glob_umin.SetLength(k + 2);

	//Vec<ZT> umin; umin.SetLength(k - j + 2);
	Vec<FT> l; l.SetLength(k - j + 2);
	Vec<FT> c; c.SetLength(k - j + 2);

	int t;

	// Set vectors
	u[j] = 1;
	l[j] = c[0] = 0;


	for(int i=j+1; i <= k; i++) {
		c[i] = l[i] = 0;
		u[i] = 0;
	}

	// Set globals
	t = j;
	glob_A = bstarnorm[0];
	glob_umin = u;

	subEnumDef<ZT, FT> new_def;
	new_def.ubar = u;
	new_def.lbar = 0;
	new_def.cbar = 0;
	new_def.tbar = t;
	L.push_back(new_def);

	while(L.size() > 0) {
		size_t len = L.size();
		subEnumDef<ZZ, RR> act_def = L[len-1];
		L.pop_back();
		cout << "Popping treeInfo" << endl;

		DagdelenSubEnumeration<ZT, FT> (B, mu, bstarnorm, A_in, act_def.tbar, k, dim, act_def);
	}
	calcVectorLength<ZT,FT>(glob_umin, B);
	return;
}




#endif /* SRC_PARALLELDAGDALEN_HPP_ */
