/*
 * old.hpp
 *
 *  Created on: 18.05.2018
 *      Author: tuxed
 */

#ifndef SRC_OLD_HPP_
#define SRC_OLD_HPP_


	/*int d = B.NumRows();
	vec_RR A; A.SetLength(d);
	for(int i=0; i<d; i++) A[i] = bstar[0] * 1.05
	// min(1.0, 1.05*(d - i) / (d - 1));
	//GamaEnumeration(mu, bstar, A);

	// Do Enmueration

	vec_ZZ u; u.SetLength(d);
	vec_ZZ u_min; u_min.SetLength(d);
	vec_RR l; l.SetLength(d + 1);
	vec_RR c; c.SetLength(d);

	// Set initial values
	for(int i=0; i<d; i++) {
		u_min[i] = 0;
		u[i] = 0;
		l[i] = 0;
		c[i] = 0;
	}
	l[d] = 0;
	u[0] = 1;
	u_min[0] = 1;

	int t = 0;


	while(t < d) {
		RR tmp_mu_sum = muProd<ZZ, RR>(u, mu, t, d-1); // Das ist c_i
		l[t] = (conv<RR>(u[t]) + tmp_mu_sum) * (conv<RR>(u[t]) + tmp_mu_sum) * bstar[t]; // Hier kommt A in erster Iteration raus
		RR lsum = vecSum(l, t, d); // Das akkumuliert der andere Algorithmus und speichert es, dann keine Wiederberechnung

		//cout << tmp_mu_sum << endl;

		RR templength = vecSum(l, 0, d);

		//cout << t << ": " << u << "("<< templength
			//	<< "," << A[0] << ")" << endl;

		if(t == 0 && templength <= A[t]) {
			u_min = u;
			u[0] += 1;
			for(int i=0; i<d; i++) A[i] = templength
					// min(1.0, 1.05*(d - i) / (d - 1));
		}

		if (t != 0 && lsum <= A[t]) {
			t--;
			RR rootval = sqrt((A[t] - vecSum(l,t+1,d)) / bstar[t]) ;
			u[t] = conv<ZZ>(ceil(-tmp_mu_sum - rootval));
		}

		if(lsum > A[t]) {
			t += 1;
			if(t >= d)
				break;
			u[t] += 1;
		}
	}

	cout << A[0] << endl;


	calcVectorLength<ZZ,RR>(u_min, B);*/


#endif /* SRC_OLD_HPP_ */
