/*
 * GramSchmidtObject.hpp
 *
 *  Created on: 28.08.2018
 *      Author: mburger
 */

#ifndef SRC_GRAMSCHMIDTOBJECT_HPP_
#define SRC_GRAMSCHMIDTOBJECT_HPP_

#define INPUT_SQUARED 0x40
#define INPUT_NONSQUARED 0x20

#include <vector>
#include <iostream>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>
#include <NTL/RR.h>

using namespace std;
using namespace NTL;

template<class T> class GramSchmidtObject {
public:
	std::vector<std::vector<T> > approxelement;
    std::vector<T> c;       //c[1..n] = |b*i|
    std::vector<T> cd;      //cd[1..n] = |b*i|^2
    std::vector<T> lena;      //len[1..n] = |b_i|^2
    std::vector<std::vector<T> >mu;      //mu[1..n][1..n]

    vec_RR cbstar;
    mat_RR mus;

    int GScomputed;

    GramSchmidtObject() {
        clear();
        GScomputed=0;
    }

    GramSchmidtObject(size_t s) {
        c.resize(s);
        cd.resize(s);
        lena.resize(s);
        cbstar.SetLength(s);

        mu.resize(s);
        mus.SetDims(s, s);

        mu.resize(s);
        for (size_t i=0;i<s;i++) {
        	mu[i].resize(s);
        }

        approxelement.resize(s);
        for (size_t i=0;i<s;i++) {
        	approxelement[i].resize(s);
        }

        GScomputed=0;
    }

    void computeGramSchmidtline(const mat_ZZ& L,std::vector<std::vector<T> >& approxelement,int index,int ioffset=1,int joffset=1) {
    	/*ZZ ip;
        int len = L[index-1].length();
        approxelement[index].resize(len+1);

        for (int k = 0;k < len; k++) {
            if (abs(L[index-1][k]) <= 1e+7) {
                approxelement[index][k+1] = to_double(L[index-1][k]);   //NTL_ZZ to T
            } else {
                conv(approxelement[index][k+1] ,L[index-1][k]);
            }

            if(k==1) c[0] = InnerProduct(B[0], B[0]);
        }

    	return;*/

        //assume GS basis for L[1..index-1] is computed
        //and compute GS for L[index]
        ZZ ip;
        int len = L[index-1].length();
        approxelement[index].resize(len+1);

        for (int i=0;i<len;i++) {
                if (abs(L[index-1][i]) <= 1e+7) {
                    approxelement[index][i+1] = to_double(L[index-1][i]);   //NTL_ZZ to T
                } else {
                    conv(approxelement[index][i+1] ,L[index-1][i]);
                }
            }

            T ipapprox;
            for (int j=1;j<=index;j++) {

                ipapprox = 0;
                int size = approxelement[index].size();

                // Inner product B[k] * B[j]
                for (int i=1;i<size;i++) {
                    ipapprox += approxelement[index][i] * approxelement[j][i];
                }

                mu[index][j] = ipapprox;

                if (j==index) {
                	lena[index] = mu[index][j];
                }

                for (int k=1;k<j;k++) {
                    mu[index][j] -= mu[j][k] * mu[index][k] * cd[k];
                }

                if (j<index) {
                    if (mu[j][j] > 0) {
                        mu[index][j] /=  mu[j][j];

                    } else {
                        mu[index][j] =0;
                    }
                }
            }

            cd[index] = mu[index][index];
            if (cd[index] > 0) {
                c[index] = sqrt(cd[index]);
            } else {
                c[index] = 0;
        }
    }

    void clear() {
        c.resize(0);
        cd.resize(0);
        lena.resize(0);
        mu.resize(0);
        approxelement.resize(0);

        GScomputed=0;
    }

    T minbi() {
        if (c.size()<=1) return 0;
        T ret = c[1];
        for (int i=1;i<c.size();i++) {
            ret = min(ret,c[i]);
        }
        return ret;
    }

    T maxbi() {
        if (c.size()<=1) return 0;
        T ret = c[1];
        for (int i=1;i<c.size();i++) {
            ret = max(ret,c[i]);
        }
        return ret;
    }

    void displaymu() {
		for (size_t i=1;i<mu.size();i++) {
			for (size_t j=1;j<mu[i].size();j++) {
				cout << mu[i][j] << " ";
			}
			cout << endl;
		}
	}

	void displaymuline(int i) {
		for (int j=1;j<=i;j++) {
			cout << mu[i][j] << " ";
		}
		cout << endl;
	}

	void displaymu(int i,int j) {
		for (int k=i;k<=j;k++) displaymuline(k);
	}

	void displayc(int istart=1,int iend=-1) {

		cout << "[";
		if (iend==-1) iend = c.size()-1;

		for (int i=istart;i<=iend;i++) {
			cout << c[i] << " ";
		}
		cout << "]" << endl;
	}

	void updateGSBasis(const mat_ZZ& L,int istart,int iend) {

		ComputeGS(L, mus, cbstar);
		//cout << "Updating GS" << endl;

		int dim = L.NumRows()-1;


		if (istart < 1) istart = 1;
		if (iend > dim) iend = dim;
		for (int i=istart;i<=iend;i++) {
			computeGramSchmidtline(L, approxelement,i);
		}
		GScomputed = iend;
	}

	void updateGSBasis(const mat_ZZ& L) {
		int n = L.NumRows() - 1;
		GScomputed = 0;
		updateGSBasis(L, 1, n);
	}

	T LatticeGH(std::vector<T>& c,char opt=INPUT_SQUARED) {
	    int i;
	    T det = T(0);
	    size_t dim = c.size()-1;
	    for (size_t i=1;i<=dim;i++) {
	        det += log(c[i]);
	    }
	    if (opt==INPUT_SQUARED) det *= 0.5;
	    det = exp(det/dim);

	    //cout << det << endl;
	    return det; // * bkzconstants::ghconstant(dim);
	}

    T LatticeGH(std::vector<T>& c,int istart,int iend,char opt=INPUT_NONSQUARED) {

    	int len = iend - istart + 1;
    	T det = T(0);

    	for(int i=istart; i <= iend; i++) {
    		det += log(cbstar[i-1]);
    	}

    	det = det / T(len);
    	det = exp(det);

    	//cout << "Det: " << det << " /";

    	//det = pow(det, 1.0 / RR(len));

    	//cout << " pow(Det): " << det << " /";

    	double gam_arg = (double)len/2.0 + 1.0;
    	gam_arg = lgamma(gam_arg);
    	T gam = exp(T(gam_arg) * T( 2.0 / double(len)));
    	gam = gam / M_PI;

    	//cout << " /Gam: " << gam << " /";

    	T frac = 1.05 * 1.05 * gam * det;
    	return frac;
    	//cout << "New Gauss: " << sqrt(frac) << endl;


        /*T ret = T(0);
        for (int i=istart;i<=iend;i++) {
            ret += log(c[i]);
        }
        ret /= (iend-istart+1);
        if (opt==INPUT_SQUARED) ret *= 0.5;
        ret = exp(ret);
        cout << ret << endl;
        return ret * 1.05;*/
    }
};

template<typename T> bool operator <(GramSchmidtObject<T>& gs1,GramSchmidtObject<T>& gs2) {
    int n = min(gs1.c.size(),gs2.c.size());
    for (int i=1;i<n;i++) {
        if (gs1.c[i] < gs2.c[i]) return true;
        if (gs1.c[i] > gs2.c[i]) return false;
    }
    return false;
}



#endif /* SRC_GRAMSCHMIDTOBJECT_HPP_ */
