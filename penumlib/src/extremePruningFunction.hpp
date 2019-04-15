#ifndef EXTREMEPRUNINGFUNCTION_HPP_
#define EXTREMEPRUNINGFUNCTION_HPP_

//#include "Utils.hpp"

#include <NTL/RR.h>
#include <NTL/LLL.h>
#include <vector>
#include <random>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <omp.h>

#include "EvenSimplex.hpp"

using boost::multiprecision::cpp_dec_float_50;


//using namespace std;
//using namespace NTL;
//using namespace boost;
using namespace std;

//NTL_CLIENT

// New trial from Matlab
template <class T>
T volBallMatlab(T R, int n) {
	T numer = pow(T(M_PI),T(n)/2.0);
	T denom = exp(boost::math::lgamma<T>(T(n)/2.0 + 1.0));
	T fac1 = pow(R, T(n));

	T res = fac1 * (numer/denom);
	return res;
}

/**
 * WORKS
 * Tested for Unitballs: Works
 */
template <class T, class FT2>
T volBall(T R, int n) {
	T fac1 = pow(T(R), T(n));

	T top = pow(T(M_PI), T(n)/T(2.0));
	T bottom = exp(boost::math::lgamma<T>(T(n)/2.0 + 1.0));

	T res = fac1 * (top/bottom);
	return res;
}

/**
 * WORKS
 * Used in randSphere()
 */
template<class T>
void sampleVectorX(T* X, int m, int n) {
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.0,1.0);

	int count = 0;
	while (count < m*n) {
		double number = distribution(generator);
			X[count] = T(number);
			count++;
	}
	return;
}

/**
 * WORKS
 *  This function returns an m by n array, X, in which
 each of the m rows has the n Cartesian coordinates
 of a random point uniformly-distributed over the
 interior of an n-dimensional hypersphere with
 radius r and center at the origin.
 */
template <class T>
void randSphere(T* X, int m, int n, T r) {
	T * Xrep = new T[m];
	T * lX = new T[m];

	sampleVectorX<T>(&X[0], m, n);

	// Calculate squaresum per
	for(int i=0; i<m; i++) {
		T sum = T(0.0);
		for(int j=0; j<n; j++) {
			sum += X[i*n + j] * X[i*n + j];
		}
		lX[i] = sum;
	}

	//int seed = std::chrono::system_clock::now().time_since_epoch().count();
	//std::default_random_engine generator (seed);
	//std::uniform_real_distribution<double> distribution(0,1.0);
	// Generate the matrix Xrep
//#pragma omp parallel for
	for(int i=0; i<m;i++) {
		T s2 = T(lX[i]);

	    //double num = distribution(generator);
	    //T gammaval = T(num);
		T gammaval = boost::math::gamma_p<T, T>(T(n)/2.0, s2/2);
		T powval = pow(gammaval, 1.0/T(n));
		T finalval = r*(powval/sqrt(lX[i]));

		Xrep[i] = finalval;

		// Scale each entry with
		for(int j=0; j<n; j++) {
			X[i*n + j] *= finalval;
		}
	}

	delete[] Xrep;
	delete[] lX;

	// cout all points sampled
	/*cout << "Res: " << endl;
	for(int i=0; i<m; i++) {
		for(int j=0; j<n; j++) {
			cout << X[i*n + j] << " ";
		}
		cout << endl;
	}*/
}



/*
 * Calculates the Volume of a cylinder intersection (R_0, ... , R_l) with Monte Carlo Sampling
 */
template <class FT1, class FT2>
FT1 volCylinderIntersection(const vector<FT2>& Ri, const FT2 radius, const int dim,
		vector<FT1>& probs, bool fillvec=true) {
	const long long m = 9*dim*dim; // Number of samples TODO: Needs to be choosen from size
	const int k = dim;
	FT1 prob;

	if(fillvec)
	{
		FT2 * X = new FT2[m * k];
		randSphere<FT2>(X, m, k, FT2(1.0));
		long long hit = 0;

		// Run over all samples
	//#pragma omp parallel for reduction(+:hit)
		for(int c=0; c < m; c++) {
			bool hitdec=true;

			// V.A. j in [1,k]
			for(int j=1; j <= k; j++) {

				// Sum from 1 to j of u_i^2
				FT2 ui2 = FT2(0.0);
				for(int i = 0; i < j; i++) {
					ui2 += X[c*k + i] * X[c*k + i];

					if(ui2 > ((Ri[j-1])*(Ri[j-1])) / ( (Ri[k-1])*(Ri[k-1]))) {
						hitdec=false;
					}
				}
			}
			if(hitdec)
				hit++;
		}
		delete[] X;
		prob = (FT1)hit/(FT1)m;
		probs[k] = prob;
	}

	else
		prob=probs[k];


	//cout << "Prob (Carlo):" << prob << endl;

	// Now caluclate the volume of the k-dimensional ball with radius radius V_k(radius)
	FT1 volB = volBall<FT1, FT2>(Ri[k-1], k);
	FT1 volCyl = volB * FT1(prob);
	return volCyl;
}

template <class FT1, class FT2>
FT1 VolumeSimplex(const vector<FT2>& Ri, const FT2 bigR, const int k, const int dim,
		vector<FT1>& probs, bool fillvec=true) {
	FT1 res = FT1(0);

	// Prepare normalized vector of prun func
	vector<FT1> vec;
	vec.push_back(0.0);

	//FT1 big = Ri[dim-1];
	for(int i=0; i<k; i+=2) {
		//FT1 tmp = Ri[i];
		//vec.push_back(tmp/big);
		vec.push_back(Ri[i]);
	}

	FT1 prob;
	if(fillvec) {
		prob = EvenSimplex::EvenSimplexVolumeWrapperMain<FT1, FT2>(vec, k/2, opt_volume_prob);
		probs[k] = prob;
		//cout << "k: " << k << " / prob: " << prob << endl;
	}
	else {
		//cout << "k: " << k << " / prob: " << probs[k] << endl;
		prob = probs[k];
	}

	//cout << "Prob: " << prob << " with k:"<< k << " and r: " << FT1(big) << endl;
	//res = prob * volBall<FT1, FT1>(sqrt(FT1(big)), k);
	res = prob * volBall<FT1, FT1>(FT1(bigR), k);
	return res;
}



template <class FT, class FT2>
FT estimateNumberOfNodesInTree (const vector<FT2>& Ri_prob, FT2 amax, const int n, const FT2* bstar, vector<FT>& probs,
		bool fillvec=true) {

	if(fillvec) {
		probs.clear();
		probs.resize(n+2);
		probs.reserve(n+2);
	}

	FT tmpsum=FT(0);

	for(int k=1; k<=n; k++) {

		FT volCyl;

		if(k%2==0 || 1==1) {
			volCyl = VolumeSimplex<FT, FT2>(Ri_prob, amax, k, k, probs, fillvec);
		}

		else
		{
            vector<FT2> Ri_root; Ri_root.reserve(n);
            for(auto it = Ri_prob.begin(); it != Ri_prob.end(); it++) {
                    Ri_root.push_back(*it * amax);
            }
            volCyl = volCylinderIntersection<FT, FT2>(Ri_root, -5, k, probs, fillvec);
		}

		FT bigprod = FT(1.0);
		for(int i=n-k; i<n; i++) {
			bigprod *= sqrt(FT(bstar[i]));
		}
		tmpsum += (volCyl/bigprod);
	}
	tmpsum *= 0.5;
	return tmpsum;
}


/**
 *  Testing to calculate the Volume within an Even Simplex
 */

template <class FT>
FT EvalAt(const vector<FT>& F,double a) {
    //Assume 0<=a<=1
	FT ret,t;
    ret = FT(0);
    t = FT(1.0);

    int i,n;
    n = F.size();
    for (i=0;i<n;i++) {
        ret += t * F[i];
        t *= FT(a);
    }
    return ret;
}

template <class FT>
vector<FT> Integ(vector<FT>& F, double low) {

	vector<FT> ret;

    int i,n;
    n = F.size();

    ret.resize(n+1); ret.reserve(n+1);
    ret[0] = 0;
    for (i=0;i<n;i++) {
        ret[i+1] = F[i] / FT(i+1);
    }
    ret[0] = -EvalAt<FT>(ret,low);

    return ret;
}


template <class FT>
FT EvenSimplexVolume(int n, double* R2, int option=0) {

    //Computing the volume of n-dimensional trancated simplex
    //{ (x1,...,xn) : \sum[i=1,...,l] x_i < R_l for l=1,...,n}
    //F[n,n] is the volume
    //pb is the probability that Pr{x <- surface of full-simplex}[x is in trancated simplex]

    int i,j;
    vector<FT>** F;
    F = new vector<FT>*[n+2];
    for (i=0;i<=n;i++) F[i] = new vector<FT>[n+2];

    F[1][0].resize(2); F[1][0].reserve(2);
    F[1][0][1] = FT(1); //F[1,0]=y

    F[1][1].resize(2); F[1][1].reserve(2);
    F[1][1][0] = FT(R2[1-1]); //F[1,0]=R_1 (stored in [0])

    for (i=2;i<=n;i++) {
        F[i][0] = Integ(F[i-1][0],0.0);
        for (j=1;j<=i;j++) {
            F[i][j] = Integ<FT>(F[i-1][j],R2[j-1]);
            F[i][j][0] += EvalAt<FT>(F[i][j-1],R2[j-1]);
            if (i==j) {
                FT v;
                v = F[i][i][0];
                for (int k=2;k<=i;k++) v *= k;
            }
        }

    }

    FT ret = F[n][n][0];
    for (i=2;i<=n;i++) ret *= i;
    for (i=0;i<=n;i++) delete [] F[i];
    delete [] F;
    //cout << "Simplex: " << ret << endl;
    return ret;
}




// Generate real number between min and max
/*template <class T>
T doubleRand(const T & min, const T & max) {
    static thread_local std::mt19937 generator;
    std::uniform_real_distribution<T> distribution(min,max);
    return distribution(generator);
}
template <class T>
void sampleRandomPoint(vector<T>& fig, const vector<T>& ref_fig, const int k) {
	for(int i = 0; i <= k; i++) {
		fig[i] = doubleRand<double>(0, ref_fig[i]);
	}
	return;
}*/

/**
 * Tested for unitballs FAILED
 */
/*template <class T>
T calc_VnR (const T R, const int n ) {
	T numerator = T(pow(M_PI, double(n)/2.0));
	T denominator = T(gamma(double(n)/2.0 + 1.0));
	T vol = T(pow(R,n)) * numerator / denominator;
	cout << "Volume is: " << vol << endl;
	return vol;
}
template <class T>
double approxCylinderVolume(const vector<T>& cyl, const vector<T>& ball, const int k) {
	vector<double> sam;
	sam.reserve(cyl.size());
	sam.resize(cyl.size());
	for(int i=0; i<1000000; i++) {
		sampleRandomPoint(sam, ball, k);
	}
	return 0.0;
}*/

#endif
