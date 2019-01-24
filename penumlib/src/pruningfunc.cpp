/*
 * pruningfunc.cpp
 *
 *  Created on: 29.06.2018
 *      Author: mburger
 */

#include "pruningfunc.hpp"
#include <iostream>
#include <math.h>

#define ALPHA 0.3

void updatePruningFuncLoc(double* prunfunc, pruning_conf& conf, const double A, const int dim,
		const int jj, const int kk) {
	if(conf.type == NO) {
		for(int i=jj; i<=kk; i++) {
			prunfunc[i] = A;
		}
	}

	else if(conf.type == LINEAR) {
		for(int i=jj; i<=kk; i++) {
			prunfunc[kk-i] = A * (i+1)/dim;

			//double val = A * min(1.0, 1.05*(dim - i) / (dim));
			//prunfunc[i] = val; // 1.05
		}
	}

	else if(conf.type == SOFT) {
		double val = A*1.05;
		for(int i=jj; i<=kk; i++) {
			prunfunc[i] = val;
			val *= 0.99;
		}
	}

	else if (conf.type == STEP) {
		for(int i=0;i<dim/2;i++) {
			prunfunc[i] = A;
		}
		for(int i=dim/2;i<dim;i++) {
			prunfunc[i] = A * ALPHA;
		}
	}

	else if(conf.type == PARLINEAR) {
		double loc_conf_param = min(conf.param, 1.0) * double(dim);
		double target = (1 + loc_conf_param) / dim;
		double step = (1.0 - target) / (double)(dim-1);
		for(int i=jj; i<=kk; i++) {
			prunfunc[/*kk-*/i] = A -
					A * (i) * step;
		}
	}

	else if(conf.type == EXTREME) {
	    double a  = -1.29185e-14;
	    double b  = 6.20066e-12;
	    double cc = -1.20165e-09;
	    double d  = 1.21218e-07;
	    double e  = -6.91288e-06;
	    double f  = 0.00022931;
	    double g  = -0.00424356;
	    double h  = 0.0400812;
	    double i  = 0.000914465;

	    double div = (double) dim;
	    double ctr = 110;
	    double y;

	    if(dim==90) {
            cout << "Special fplll pruning func 90" << endl;
            double fac[91];
            fac[90] = 1;
            fac[89] = 1;
			fac[88] = 1;
			fac[87] = 1;
			fac[86] = 1;
			fac[85] = 1;
			fac[84] = 0.999554482;
			fac[83] = 0.999554482;
			fac[82] = 0.993081058;
			fac[81] = 0.993081058;
			fac[80] = 0.992247743;
			fac[79] = 0.992247743;
			fac[78] = 0.976992703;
			fac[77] = 0.976992703;
			fac[76] = 0.959580392;
			fac[75] = 0.959580392;
			fac[74] = 0.941544905;
			fac[73] = 0.941544905;
			fac[72] = 0.920816028;
			fac[71] = 0.920816028;
			fac[70] = 0.899402344;
			fac[69] = 0.899402344;
			fac[68] = 0.869385651;
			fac[67] = 0.869385651;
			fac[66] = 0.837911108;
			fac[65] = 0.837911108;
			fac[64] = 0.803662835;
			fac[63] = 0.803662835;
			fac[62] = 0.767891426;
			fac[61] = 0.767891426;
			fac[60] = 0.731804441;
			fac[59] = 0.731804441;
			fac[58] = 0.695679801;
			fac[57] = 0.695679801;
			fac[56] = 0.659637217;
			fac[55] = 0.659637217;
			fac[54] = 0.624299232;
			fac[53] = 0.624299232;
			fac[52] = 0.59111711;
			fac[51] = 0.59111711;
			fac[50] = 0.559071864;
			fac[49] = 0.559071864;
			fac[48] = 0.528427676;
			fac[47] = 0.528427676;
			fac[46] = 0.49962314;
			fac[45] = 0.49962314;
			fac[44] = 0.471429732;
			fac[43] = 0.471429732;
			fac[42] = 0.445024764;
			fac[41] = 0.445024764;
			fac[40] = 0.420964399;
			fac[39] = 0.420964399;
			fac[38] = 0.397456438;
			fac[37] = 0.397456438;
			fac[36] = 0.373815546;
			fac[35] = 0.373815546;
			fac[34] = 0.352107282;
			fac[33] = 0.352107282;
			fac[32] = 0.333496135;
			fac[31] = 0.333496135;
			fac[30] = 0.315557915;
			fac[29] = 0.315557915;
			fac[28] = 0.300254085;
			fac[27] = 0.300254085;
			fac[26] = 0.286103739;
			fac[25] = 0.286103739;
			fac[24] = 0.27334306;
			fac[23] = 0.27334306;
			fac[22] = 0.263000691;
			fac[21] = 0.263000691;
			fac[20] = 0.257205086;
			fac[19] = 0.257205086;
			fac[18] = 0.248183482;
			fac[17] = 0.248183482;
			fac[16] = 0.235154146;
			fac[15] = 0.235154146;
			fac[14] = 0.219664612;
			fac[13] = 0.219664612;
			fac[12] = 0.202518238;
			fac[11] = 0.202518238;
			fac[10] = 0.184161696;
			fac[9] = 0.184161696;
			fac[8] = 0.164826542;
			fac[7] = 0.164826542;
			fac[6] = 0.144579344;
			fac[5] = 0.144579344;
			fac[4] = 0.12339341;
			fac[3] = 0.12339341;
			fac[2] = 0.10220174;
			fac[1] = 0.10220174;

        	for (int ii = 0; ii<dim; ii++){
                // To avoid rounding errors because in the polynom, it is not always 1.00 * A
                fac[ii] = min((conf.param) + fac[ii+1], 1.0);
                prunfunc[dim-ii-1] = A * fac[ii+1];
                //cout << fac[ii] << " ";
            }
	    }
        
        if(dim==66) {
            cout << "Special fplll pruning func 66" << endl;
            double fac[66];
            fac[65] = 1;
            fac[64] = 1;
            fac[63] = 1;
            fac[62] = 1;
            fac[61] = 0.999998743;
            fac[60] = 0.999998743;
            fac[59] = 0.999632281;
            fac[58] = 0.999632281;
            fac[57] = 0.997742402;
            fac[56] = 0.997742402;
            fac[55] = 0.983964966;
            fac[54] = 0.983964966;
            fac[53] = 0.959440876;
            fac[52] = 0.959440876;
            fac[51] = 0.927132817;
            fac[50] = 0.927132817;
            fac[49] = 0.887801533;
            fac[48] = 0.887801533;
            fac[47] = 0.844829554;
            fac[46] = 0.844829554;
            fac[45] = 0.801179895;
            fac[44] = 0.801179895;
            fac[43] = 0.759059704;
            fac[42] = 0.759059704;
            fac[41] = 0.71605206;
            fac[40] = 0.71605206;
            fac[39] = 0.671886668;
            fac[38] = 0.671886668;
            fac[37] = 0.628687347;
            fac[36] = 0.628687347;
            fac[35] = 0.586685552;
            fac[34] = 0.586685552;
            fac[33] = 0.544654333;
            fac[32] = 0.544654333;
            fac[31] = 0.506013112;
            fac[30] = 0.506013112;
            fac[29] = 0.470508894;
            fac[28] = 0.470508894;
            fac[27] = 0.438618024;
            fac[26] = 0.438618024;
            fac[25] = 0.410204712;
            fac[24] = 0.410204712;
            fac[23] = 0.38412724;
            fac[22] = 0.38412724;
            fac[21] = 0.361383161;
            fac[20] = 0.361383161;
            fac[19] = 0.340633166;
            fac[18] = 0.340633166;
            fac[17] = 0.321921202;
            fac[16] = 0.321921202;
            fac[15] = 0.305365541;
            fac[14] = 0.305365541;
            fac[13] = 0.290098303;
            fac[12] = 0.290098303;
            fac[11] = 0.274343704;
            fac[10] = 0.274343704;
            fac[9] = 0.254620279;
            fac[8] = 0.254620279;
            fac[7] = 0.230806954;
            fac[6] = 0.230806954;
            fac[5] = 0.203306703;
            fac[4] = 0.203306703;
            fac[3] = 0.171874674;
            fac[2] = 0.171874674;
            fac[1] = 0.13455723;
            fac[0] = 0.13455723;            
        
        	for (int ii = 0; ii<dim; ii++){
                // To avoid rounding errors because in the polynom, it is not always 1.00 * A
                fac[ii] = min((conf.param) + fac[ii], 1.0);
                prunfunc[dim-ii-1] = A * fac[ii];
                //cout << fac[ii] << " ";
            }
	    }
        
        else {
            for (int ii = 0; ii<dim; ii++){
	    	// To avoid rounding errors because in the polynom, it is not always 1.00 * A
	    	y = (double)ii * ctr / div;
	    	double fac = min(1.0, ( a * pow(y, 8) + b *pow(y,7) + cc * pow(y,6) + d * pow(y,5) + e * pow(y,4) + f * pow(y,3) + g * pow(y,2) + h * y+i));
	    	fac = min((conf.param) + fac, 1.0);
	    	prunfunc[dim-ii-1] = A * fac;
            //cout << fac << " ";
            }
        }
    }

	// Default is no pruning
	else {
		for(int i=0; i<dim; i++) {
			prunfunc[i] = A;
		}
	}
	prunfunc[kk+1] = prunfunc[kk];
}
