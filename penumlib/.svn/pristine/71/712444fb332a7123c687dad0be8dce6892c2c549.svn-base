/*
 * pruningfunc.hpp
 *
 *  Created on: 29.06.2018
 *      Author: mburger
 */

#ifndef SRC_PRUNINGFUNC_HPP_
#define SRC_PRUNINGFUNC_HPP_

#include "MBVec.hpp"

enum Pruning {NO, LINEAR, SOFT, STEP, PARLINEAR, EXTREME};

struct pruning_conf {
	Pruning type;
	double param;
};

struct pruning_conf2{
	int type;
	double param;
} ;

void updatePruningFuncLoc(double* prunfunc, pruning_conf& conf, const double A, const int dim,
		const int jj, const int kk);


#endif /* SRC_PRUNINGFUNC_HPP_ */
