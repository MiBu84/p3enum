#ifndef EXTREMEPRUNINGFUNCTION_HPP_
#define EXTREMEPRUNINGFUNCTION_HPP_

#include "Utils.hpp"

#include <NTL/RR.h>
#include <NTL/LLL.h>
#include <vector>
#include <random>
#include <boost/math/special_functions/gamma.hpp>


//using namespace std;
using namespace NTL;
using namespace boost;

NTL_CLIENT

template<class T>
void sampleVectorX(double* X, int m, int n) {
	std::default_random_engine generator;
	std::normal_distribution<double> distribution(0.0,1.0);

	int count = 0;
	while (count < m*n) {
		double number = distribution(generator);
			X[count] = number;
			count++;
	}
	return;
}

template <class T>
void randSphere(double* X, int m, int n, double r) {
	double * Xrep = new double[m];
	double * lX = new double[m];

	sampleVectorX<T>(&X[0], m, n);

	// Calculate squaresum per
	for(int i=0; i<m; i++) {
		double sum = 0.0;
		for(int j=0; j<n; j++) {
			sum += X[i*n + j] * X[i*n + j];
		}
		lX[i] = sum;
	}

	// Generate the matrix Xrep
	for(int i=0; i<m;i++) {
		double s2 = lX[i];

		double gammaval = math::gamma_p<double, double>(double(n)/2.0, s2/2);
		double powval = pow(gammaval, 1.0/double(n));
		double finalval = r*(powval/sqrt(lX[i]));

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

template <class T>
double volBall(double R, int n) {
	double fac1 = pow(R, n);

	double top = pow(M_PI, double(n)/2.0);
	double bottom = exp(lgamma(double(n)/2.0 + 1.0));

	return fac1 * (top/bottom);
}


template <class T>
double volCylinderIntersection(const vector<double>& Ri, const double radius) {
	const int m = 2000;
	const int n = 3;
	T * X = new T[m * n];

	randSphere<T>(X, m, n, radius);

	int hit = 0;

	for(int k=0; k<m; k++) {
		// Check condition
		bool hitdec=true;
		for(int j = 0; j < n; j++) {

			// Call Sum of u_i^2
			double ui2 = 0.0;
			for(int i = 0; i<=j; i++) {
				ui2 += X[k*n + i] * X[k*n + i];
			}
			if(ui2 >(Ri[j]*Ri[j]) / ( Ri[n-1]*Ri[n-1])  ) {
				hitdec=false;
			}
		}
		if(hitdec)
			hit++;
	}


	cout << "We had " << hit << " hits." << endl;
	delete[] X;
	return 0.0;
}

// Generate real number between min and max
template <class T>
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
}

template <class T>
T calc_VnR (const T R, const int n ) {

	T numerator = T(pow(M_PI, double(n)/2.0));
	T denominator = T(gamma(double(n)/2.0 + 1.0));
	T vol = T(pow(R,RR(n))) * numerator / denominator;

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
}

#endif
