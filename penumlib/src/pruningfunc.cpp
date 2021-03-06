/*
 * pruningfunc.cpp
 *
 *  Created on: 29.06.2018
 *      Author: mburger
 */

#include "pruningfunc.hpp"
#include "Configurator.hpp"
#include "Utils.hpp"
#include <iostream>
#include <string>
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

		//cout << "A is " << A << endl;

		if (Configurator::getInstance().ext_pruning_function_file.compare("NOT") != 0 ) {
			const vector<double>& fac = Configurator::getInstance().ext_pruning_function;

			if((int)fac.size() != dim) {
				cout << "Size mismatch between pruning and lattice: " <<  fac.size() << " / " << dim << endl;
				exit(-55);
			}

			//cout << "Read prun func." << endl;
        	for (int ii = 0; ii<dim; ii++){
                // To avoid rounding errors because in the polynom, it is not always 1.00 * A                
                prunfunc[dim-ii-1] = A * min((conf.param) + fac[ii], 1.0);
                //cout << "PrunFunc[" << dim-ii-1 << "] = " << prunfunc[dim-ii-1] << endl;
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

    	/*for (int ii = 0; ii<dim; ii++){
			cout << prunfunc[ii] << endl;
        }
		exit(-1);*/
    }

	// Default is no pruning
	else {
		for(int i=0; i<dim; i++) {
			prunfunc[i] = A;
		}
	}
	prunfunc[kk+1] = prunfunc[kk];
}
