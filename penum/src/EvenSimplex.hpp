/*
 * EvenSimplex.hpp
 *
 *  Created on: 15.02.2019
 *      Author: mburger
 */

#ifndef SRC_EVENSIMPLEX_HPP_
#define SRC_EVENSIMPLEX_HPP_

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <cmath>

#include "SharedMemory.hpp"

using namespace boost::multiprecision;

#define opt_upper 0x01
#define opt_lower 0x02

#define opt_volume 0x01
#define opt_volume_prob 0x02
#define opt_surface 0x03
#define opt_surface_prob 0x04
#define opt_gaussian_prob 0x05

typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<10> > float10;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<15> > float15;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<20> > float20;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<25> > float25;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<30> > float30;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<35> > float35;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<40> > float40;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<45> > float45;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<50> > float50;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<80> > float80;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100> > float100;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<150> > float150;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<200> > float200;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<300> > float300;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<400> > float400;
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<500> > float500;

template <typename T> class autoarray1d {
    std::vector<T> data;

    public:
    T& operator [] (int i) {
        if (i<0) {
            cerr << "autoarray1d: range error" << endl;
            exit(0);
        }

        if ((cpp_int)i>=data.size()) {
            data.resize(i+1);
        }
        return data[i];
    }
    int size() {return data.size();};
    void resize(int i) {data.resize(i);};

    void load(std::string& fname) {
        //loadstdvector(data,fname);
    }

    void save(std::string& fname) {
        //savestdvector(data,fname);
    }
};

template <typename T> class autoarray2d {
    autoarray1d<autoarray1d<T> > data;

    public:
    autoarray1d<T>& operator [] (int i) {
        if (i<0) {
            cerr << "autoarray2d: range error" << endl;
            exit(0);
        }

        if (i>=data.size()) {
            data.resize(i+1);
        }
        return data[i];
    }
    int size() {return data.size();};
    void resize(int i) {data.resize(i);};

    void load(std::string& fname) {
//
    }

    void save(std::string& fname) {
//
    }

};

namespace bkzconstants {

bool initialized = false;
autoarray1d<cpp_int> factorial_table;
autoarray2d<cpp_int> binary_coeff_table;


void initialize() {
    if (initialized == false) {
        initialized = true;
    }
}

/**
 * binary_coeff
 * Returns positive non-zero values. Is this correct?
 */
cpp_int binary_coeff(int n,int k) {

        if (k<=0) return 1;
        if (k==1) return n;
        if (k>n) return 1;
        if (n==1) return 1;

        if (k > n/2) return binary_coeff(n,n-k);

        if ((n>100000) || (k>100000)) {
            cerr << "binary_coeff error: n=" << n << " k=" << k << endl;
            exit(0);
        }

        int cs = binary_coeff_table.size();
        if (n >= cs) {
            #pragma omp critical
            {
                binary_coeff_table.resize(n+1);
            }
        }

        cs = binary_coeff_table[n].size();
        if (k >= cs) {
            #pragma omp critical
            {
                binary_coeff_table[n].resize(n/2+2);
            }
        }

        if (binary_coeff_table[n][k]==0) {
            #pragma omp critical
            {
                binary_coeff_table[n][0] = 1;
                binary_coeff_table[n][1] = n;
                for (int kk=2;kk<=n/2;kk++) {
                    binary_coeff_table[n][kk] = binary_coeff_table[n][kk-1] * (n-kk+1)/kk;
                    //cout << "n=" << n << " k=" << kk << " " << binary_coeff_table[n][kk-1]  << " " << binary_coeff_table[n][kk] << endl;
                }
            }
        }

        return binary_coeff_table[n][k];
    }

/**
 * Factorial
 * works good
 */
cpp_int factorial(int k) {

    initialize();
    if (k<=0) return 0;
    if (k==1) return 1;

    if (k>100000) {
        cerr << "fractional: k error " << k << endl;
        exit(0);
    }

    if (factorial_table.size()<3) {
        factorial_table.resize(3);
        factorial_table[1] = 1;
        factorial_table[2] = 2;
    }
    int cs = factorial_table.size();
    if (k >= cs) {
        #pragma omp critical
        {
            factorial_table.resize(k+1);
            for (int i=cs;i<=k;i++) {
            	factorial_table[i] = factorial_table[i-1] * i ;
            }
        }
    }
    return factorial_table[k];
}
} // End namespace bkzconstants


namespace EvenSimplex {

	template <typename T> int GetESprecision(int dim, T uprob) {
		//uprob is guaranteed probability bound
		int prec = 1 + dim / 2.1;   //for full-enum
		if (uprob <= 0.1) {
			prec = std::min(prec,(int)(2 + dim/3.6));
		}
		if (prec < 4) return 4;
		return prec;
	}


/**
 *
 */

	template <typename DFLOAT,typename CFLOAT> DFLOAT EvenSimplexVolumeWrapperMain(std::vector<DFLOAT>& rd,int dim,int option,int parallel=1,int usedprecision=-1,DFLOAT uprob=1,DFLOAT* diag=0) {

	    //Computing the volume of n-dimensional trancated simplex
	    //{ (x1,...,xn) : \sum[i=1,...,l] x_i < R_l for l=1,...,n}
	    //F[n,n] is the volume
	    //pb is the probability that Pr{x <- surface of full-simplex}[x is in trancated simplex]

	    if (dim==1) return rd[1];

	    //SetPrecision
	    int prec = 0;
	    int precisionback = 0;

	    bool mpfr_flag = false;
	    if ((typeid(CFLOAT) != typeid(double)) && (typeid(CFLOAT) != typeid(long double)))  {
	        if (usedprecision==-1) {
	            prec = EvenSimplex::GetESprecision(dim,uprob);
	        } else {
	            prec = usedprecision;
	        }
	        precisionback = mpfr_float::default_precision() ;
	        mpfr_float::default_precision(prec);
	        mpfr_flag = true;
	    }

	    CFLOAT** F = (CFLOAT**)shared_memory::allocate2<CFLOAT>(14002,dim+1,dim+1);
	    //CFLOAT* bin = (CFLOAT*)shared_memory::allocate1<CFLOAT>(3332,dim+1);
	    CFLOAT* CR = (CFLOAT*)shared_memory::allocate1<CFLOAT>(3333,dim+1);
	    for (int i=1;i<=dim;i++) CR[i] = boost::lexical_cast<CFLOAT>(rd[i]);
	    for (int i=0;i<dim;i++) F[0][i] = 0;

	    if (option==opt_volume_prob) {
	        F[0][dim] = 1.0;   //computes F[d,d] * d!, probability
	    } else
	    if (option==opt_surface_prob) {
	        F[0][dim] = 1.0/dim;   //computes F[d,d] * d!, probability
	    }
	    if (option==opt_volume) {
	        F[0][dim] =  (CFLOAT)1.0/(CFLOAT)bkzconstants::factorial(dim);
	    } else
	    if (option==opt_surface) {
	        F[0][dim] =  (CFLOAT)1.0/(CFLOAT)bkzconstants::factorial(dim);
	    }

	    for (int j=1;j<=dim;j++) {
	        CFLOAT ht = F[j-1][dim-(j-1)];
	        CFLOAT pw = 1;
	        for (int k=dim-j+1;k>=0;k--) {
	            F[j][k] = F[j-1][k] - ht * pw * (CFLOAT)bkzconstants::binary_coeff(dim-j+1,k);
	            //cout << F[j][k] << endl;
	            if (k>0) pw *= (-CR[j]);
	        }
	    }

	    for (int i=1;i<dim;i++) {
	        F[i][0] = F[i][dim-i] * (CFLOAT)bkzconstants::factorial(dim-i);
	    }
	    if (diag!=0) {
	        for (int i=1;i<=dim;i++) diag[i] = boost::lexical_cast<DFLOAT>(F[i][0]);
	    }
	    DFLOAT ret;

	    if ((option==opt_volume) || (option==opt_volume_prob)) {
	        ret = boost::lexical_cast<DFLOAT>(F[dim][0]);
	    }
	    if ((option==opt_surface) || (option==opt_surface_prob)) {
	        ret = boost::lexical_cast<DFLOAT>(F[dim-1][0]);
	    }

	    //if CFLOAT==mpfr_float
	    if (mpfr_flag == true) {
	        mpfr_float::default_precision(precisionback);
	    }
	    return ret;
	}




/**
 *
 */
	template <typename DFLOAT> DFLOAT EvenSimplexVolumeWrapper(std::vector<DFLOAT> rd,int dim,int option,int parallel,DFLOAT uprob=1,DFLOAT* diag=0) {
	    if (uprob > 0.1) {
	        if (dim<=35) {
	            return EvenSimplexVolumeWrapperMain<DFLOAT,double>(rd,dim,option,parallel,-1,uprob,diag);
	        } else
	        if (dim<=40) {
	            return EvenSimplexVolumeWrapperMain<DFLOAT,long double>(rd,dim,option,parallel,-1,uprob,diag);
	        } else
	        if (dim<=85) {
	            return EvenSimplexVolumeWrapperMain<DFLOAT,float10>(rd,dim,option,parallel,-1,uprob,diag);
	        } else
	        if (dim<=115) {
	            return EvenSimplexVolumeWrapperMain<DFLOAT,float30>(rd,dim,option,parallel,-1,uprob,diag);
	        } else
	        if (dim<=165) {
	            return EvenSimplexVolumeWrapperMain<DFLOAT,float50>(rd,dim,option,parallel,-1,uprob,diag);
	        } else
	        if (dim<=215) {
	            return EvenSimplexVolumeWrapperMain<DFLOAT,float80>(rd,dim,option,parallel,-1,uprob,diag);
	        } else
	        if (dim<=265) {
	            return EvenSimplexVolumeWrapperMain<DFLOAT,float100>(rd,dim,option,parallel,-1,uprob,diag);
	        } else
	        if (dim<=360) {
	            return EvenSimplexVolumeWrapperMain<DFLOAT,float150>(rd,dim,option,parallel,-1,uprob,diag);
	        } else
	        if (dim<=460) {
	            return EvenSimplexVolumeWrapperMain<DFLOAT,float200>(rd,dim,option,parallel,-1,uprob,diag);
	        } else
	        if (dim<=700) {
	            return EvenSimplexVolumeWrapperMain<DFLOAT,float300>(rd,dim,option,parallel,-1,uprob,diag);
	        } else {

	            int usedprecision = EvenSimplex::GetESprecision(dim,uprob);
	            bkzconstants::factorial(dim+1);
	            DFLOAT ret;
	            #pragma omp critical
	            {
	                ret = EvenSimplexVolumeWrapperMain<DFLOAT,mpfr_float>(rd,dim,option,parallel,usedprecision,uprob,diag);
	            }
	            return ret;
	            }
	    } else {
	        if (dim<=54) {
	            return EvenSimplexVolumeWrapperMain<DFLOAT,double>(rd,dim,option,parallel,-1,uprob,diag);
	        } else
	        if (dim<=70) {
	            return EvenSimplexVolumeWrapperMain<DFLOAT,long double>(rd,dim,option,parallel,-1,uprob,diag);
	        } else
	        if (dim<=100) {
	            return EvenSimplexVolumeWrapperMain<DFLOAT,float10>(rd,dim,option,parallel,-1,uprob,diag);
	        } else
	        if (dim<=160) {
	            return EvenSimplexVolumeWrapperMain<DFLOAT,float30>(rd,dim,option,parallel,-1,uprob,diag);
	        } else
	        if (dim<=260) {
	            return EvenSimplexVolumeWrapperMain<DFLOAT,float50>(rd,dim,option,parallel,-1,uprob,diag);
	        } else
	        if (dim<=440) {
	            return EvenSimplexVolumeWrapperMain<DFLOAT,float100>(rd,dim,option,parallel,-1,uprob,diag);
	        } else
	        if (dim<=620) {
	            return EvenSimplexVolumeWrapperMain<DFLOAT,float150>(rd,dim,option,parallel,-1,uprob,diag);
	        } else
	        if (dim<=800) {
	            return EvenSimplexVolumeWrapperMain<DFLOAT,float200>(rd,dim,option,parallel,-1,uprob,diag);
	        } else
	        if (dim<=1160) {
	            return EvenSimplexVolumeWrapperMain<DFLOAT,float300>(rd,dim,option,parallel,-1,uprob,diag);
	        } else {
	            int usedprecision = EvenSimplex::GetESprecision(dim,uprob);
	            DFLOAT ret;
	            bkzconstants::factorial(dim+1);
	            #pragma omp critical
	            {
	                ret =  EvenSimplexVolumeWrapperMain<DFLOAT,mpfr_float>(rd,dim,option,parallel,usedprecision,uprob,diag);
	            }
	            return ret;
	        }
	    }

	}
}


#endif /* SRC_EVENSIMPLEX_HPP_ */
