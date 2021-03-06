//============================================================================
// Name        : HelloWord.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

//#include "parallelDagdalen.hpp"
#include "pEnumeratorDouble.hpp"
#include "BurgerEnum.hpp"
#include "Configurator.hpp"
#include "Utils.hpp"
#include "GramSchmidtObject.hpp"

#include <omp.h>
#include <NTL/mat_RR.h>
#include <cmath>
#include <string.h>

#include "extremePruningFunction.hpp"

//#include <Eigen/Dense>
using namespace std;
//using namespace Eigen;

using namespace std;
using namespace NTL;


/*template <class Matrix>
Matrix gramschmidt( const Matrix & A ) {
    Matrix Q = A;
    // First vector just gets normalized
    Q.col(0).normalize();
    for(unsigned int j = 1; j < A.cols(); ++j) {
        // Replace inner loop over each previous vector in Q with fast matrix-vector multiplication
        Q.col(j) -= Q.leftCols(j) * (Q.leftCols(j).transpose() * A.col(j));
        // Normalize vector if possible (othw. means colums of A almsost lin. dep.
        if( Q.col(j).norm() <= 10e-14 * A.col(j).norm() ) {
            std::cerr << "Gram-Schmidt failed because A has lin. dep columns. Bye." << std::endl;
            break;
        } else {
            Q.col(j).normalize();
        }
    }
    return Q;
}*/

template <class ZT>
vec_RR vec_normalize(NTL::Vec<ZT>& vec) {
	RR sum = RR(0);

	for(int i=0; i< vec.length(); i++) {
		sum += vec[i] * vec[i];
	}

	RR len = NTL::SqrRoot(sum);

	for(int i=0; i< vec.length(); i++) {
		vec[i] = vec[i] / len;
	}
	return vec;
}


/*mat_RR GramSchmidtZZ( const mat_ZZ A ) {

	mat_RR Arr; Arr.SetDims(A.NumRows(), A.NumCols());
	mat_RR Q; Q.SetDims(A.NumRows(), A.NumCols());

	for(unsigned int r = 0; r < A.NumRows(); ++r)
		for(unsigned int c = 0; c < A.NumCols(); ++c)
			Arr[r][c] = RR(A[r][c]);

    // First vector just gets normalized
	//for(unsigned int n = 0; n < A.NumCols(); ++n) {
		//Q[0][n] = RR(A[0][n]);
	//}

	Q[0] = Arr[0];
    Q[0] = vec_normalize(Q[0]);

    for(unsigned int n = 1; n < A.NumRows(); ++n) {
    	vec_RR new_col; new_col.SetLength(A.NumCols());
    	new_col = A[n];
    	for(unsigned int i=0; i <= n-1; i++ ) {
    		new_col -= (Q[i] * A[n]) / (Q[i] * Q[i]) * Q[i];
    	}

    }


   for(unsigned int j = 1; j < A.cols(); ++j) {
        // Replace inner loop over each previous vector in Q with fast matrix-vector multiplication
        Q.col(j) -= Q.leftCols(j) * (Q.leftCols(j).transpose() * A.col(j));


        // Normalize vector if possible (othw. means colums of A almsost lin. dep.
        if( Q.col(j).norm() <= 10e-14 * A.col(j).norm() ) {
            std::cerr << "Gram-Schmidt failed because A has lin. dep columns. Bye." << std::endl;
            break;
        } else {
            Q.col(j).normalize();
        }
    }
    return Q;
}*/


/*template<class ZT, class FT> FT muProd (const NTL::Vec<ZT>& x, const NTL::Mat<FT>& mu, const int t, const int s) {
	FT sum; sum=0;

	for(int i = t + 1; i <= s; i++) {
		FT convx = conv<FT>(x[i]);
		sum += convx * conv<FT>(mu[i][t]);
	}

	return sum;
}

// Calculate length
template<class ZT, class FT> FT calcVectorLength(const NTL::Vec<ZT>& coeffs, const NTL::Mat<ZT>& B) {

	int d = B.NumCols();
	cout << "coefficients: [ ";
	for(int i=0; i < d; i++) {
		cout << coeffs[i] << " ";
	}
	cout << "]" << endl;

	ZT len; len = 0;
	NTL::Vec<ZT> shortvec; shortvec.SetLength(d);

	for(int c = 0; c < d; c++) {
		shortvec[c] = 0;
	}

	cout << "Shortest Vector: [";
	for(int c = 0; c < d; c++) {
		for(int r = 0; r < d; r++) {
			shortvec[c] += coeffs[r] * conv<ZT>(B[r][c]);
		}
		cout << shortvec[c] << " ";

	}
	cout << "]" << endl;

	for(int i=0; i < d; i++) {
		len += shortvec[i] * shortvec[i];
	}
	cout << "Len: " << sqrt(conv<FT>(len)) << endl;
	cout << endl;
	return conv<FT>(len);
}

template <class ZT, class FT> void basisQualityAnalysis(const NTL::Mat<ZT>& B) {
	NTL::Mat<FT> mu;
	NTL::Vec<FT> bstar;
	ComputeGS(B, mu, bstar);

	cout << "Len of ||b_o|| =  " << bstar[0] << endl;;

	FT av_slope = FT(0.0);
	FT min_slope = FT(100000000000.0);
	FT max_slope = FT(0.0);

	std::ofstream outfile;

    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    //oss << std::put_time(&tm, "%d-%m-%Y_%H-%M-%S");
    auto str = oss.str();
    std::string filename = "Bdiag.txt" + str;

	outfile.open(filename, std::ios_base::out);
	cout << "Writing to file " << filename << endl;

	for(int i=0; i<B.NumCols(); i++) {
		outfile << i  << "\t" << log10(bstar[i]) << endl;
	}

	outfile.close();

	for(int i=0; i<B.NumCols()-1; i++) {
		RR slope = 1.0 / (log10(bstar[i+1]) - log10(bstar[i]));
		av_slope+=slope;
		if(slope>max_slope)
			max_slope=slope;
		if(slope<min_slope)
			min_slope=slope;
	}

	av_slope/=B.NumCols()-1;
	cout << "Av. slope: " << av_slope
			<< " / min. slope " << min_slope
			<< " / max. slope " << max_slope << endl;

	return;
}



RR vecSum(const vec_RR& l, const int start, const int d) {
	RR sum; sum = 0;

	for(int j = start; j < d; j++) {
		sum += l[j];
	}
	return sum;
}


template <typename VEC> void printVector (const VEC& vec) {
	for(auto it=vec.begin(); it!=vec.end(); it++) {
		cout << *it << " ";
	}
	cout << endl;
}





template<class ZT, class FT> FT BurgerEnumeration(const NTL::Mat<ZT>& B, const NTL::Mat<FT>& mu, const NTL::Vec<FT>& bstarnorm,
		NTL::Vec<ZT>& u, const FT& A_in, int j, int k, int dim, int mode=0) {

	// When subtrees are enumerated, then all possible sign-combinations must be checked
	//bool fix_sign = !(k < dim - 1);
	bool fix_sign = !(mode==2);

	//cout << "Sign check: " << fix_sign << endl;

	NTL::Vec<ZT> umin; umin.SetLength(dim + 1);//(k - j + 2);
	umin = u;
	NTL::Vec<ZT> v; v.SetLength(dim + 1);//(k - j + 2);
	NTL::Vec<FT> l; l.SetLength(dim + 1);//(k - j + 2);
	NTL::Vec<FT> c; c.SetLength(dim + 1);//(k - j + 2);

	// To decide the zigzag pattern
	NTL::Vec<ZT> Delta; Delta.SetLength(dim + 1);//(k - j + 2);
	NTL::Vec<ZT> delta; delta.SetLength(dim + 1);//(k - j + 2);

	for(int i=0; i <= B.NumCols(); i++) {
		c[i] = l[i] = 0;
		Delta[i] = v[i] = 0;
		delta[i] = 1;
	}

	int s, t; // s required for zigzagpattern
	s = t = j;

	// initial calculation of c and l if some parts of the tree
	// are already processed
	for(int i = j; i < dim; i++) {
		c[i] = muProd<ZT, FT>(u, mu, i, dim-1);
	}

	for(int i=dim-1; i>=0; i--) {
		l[i] = l[ i + 1 ] + ( conv<FT>(u[i]) + c[i] ) * ( conv<FT>(u[i]) + c[i] ) * bstarnorm[i];
	}

	// It is impossible that length is zero for non-zero vector
	// To prevent finding zero vector, then set [1, 0, 0, ... 0]
	if (l[0] == 0 && j == 0)
		u[0] = 1;
	umin = u;

	// Set scalars
	FT A = A_in;

	while (t <= k) {
		l[t] = l[ t + 1 ] + ( conv<FT>(u[t]) + c[t] ) * ( conv<FT>(u[t]) + c[t] ) * bstarnorm[t];
		//if (t==0)
		//	cout << "u: " << u << " / l: " << l << " / c: " << c << "@ t=" << t << " with A=" << l[t] << endl;
		if(l[t] < A) {
			if(t > j) {
				t = t - 1;
				// Hier aendert sich jetzt etwas fuer alles darunter liegende!
				c[t] = muProd<ZT, FT>(u, mu, t, dim-1);
				u[t] = v[t] = conv<ZT>(ceil(-c[t] - 0.5));

				// For ZigZag
				Delta[t] = 0;
				if(conv<FT>(u[t]) > -c[t]) {
					delta[t] = - 1;
				}
				else {
					delta[t] = 1;
				}
			}

			else {
				if(mode <= 1) {
					A = l[t];
					umin = u;
				}
				else if (mode == 2)
					glob_candidates_queue.push(u);

				// Just find one short vector
				if(mode==1) {
					return A;
				}

				// Avoid visiting short vector twice
				// When bound is not updated
				t = t + 1;
				s = max(s,t);

				if(!fix_sign)
					Delta[t] = -Delta[t];
				else if (t<s)
					Delta[t] = -Delta[t];

				if(Delta[t] * delta[t] >= 0) Delta[t] += delta[t];

				u[t] = v[t] + Delta[t];
			}
		}


		else {
			t = t + 1;
			s = max(s,t);

			if(!fix_sign)
				Delta[t] = -Delta[t];
			else if (t<s)
				Delta[t] = -Delta[t];

			if(Delta[t] * delta[t] >= 0) Delta[t] += delta[t];

			u[t] = v[t] + Delta[t];
		}

	}
	u = umin;
	//cout << "A: " << A << endl;
	return A;

	s = max(s,t); //

	if(t < s) {
		if(Delta[t] == 0) Delta[t]++;
		else if(Delta[t] < 0) Delta[t] = abs(Delta[t]) + 1;
		else if(Delta[t] > 0) Delta[t] = - Delta[t];
	}
	else {
		Delta[t]++;
	}
	u[t] = v[t] + Delta[t];
}

template <typename FT, typename ZT> void GamaEnumeration(const NTL::Mat<ZT>& B, const NTL::Mat<FT>& mu, const NTL::Vec<FT>& bstarnorm, const FT& A) {
	int dim = bstarnorm.length();

	NTL::Mat<FT> sigma; sigma.SetDims(dim+1,dim+1);
	NTL::Vec<FT> r; r.SetLength(dim+1);
	NTL::Vec<FT> rho; rho.SetLength(dim+1); // partial norm
	NTL::Vec<FT> v; v.SetLength(dim+1); // das ist u
	NTL::Vec<FT> v_min; v_min.SetLength(dim+1); // das ist u
	NTL::Vec<FT> c; c.SetLength(dim+1); // centers
	NTL::Vec<FT> w; w.SetLength(dim+1); // jumps

	// All vectors

	for(int i=0; i<=dim; i++) {
		r[i] = dim - 1;
		rho[i] = 0;
		v[i] = 0;
		c[i] = 0;
		w[i] = 0;
	}

	r[0] = 0;
	v[0] = 1;

	// All matrices
	for(int i=0; i<=dim; i++)
		for(int j=0; j<=dim; j++)
			sigma[i][j] = 0;

	// All scalars
	int last_nonzero = 0;
	int k = 0;
	FT R = A;

	while (true) {
		rho[k] = rho[k+1] + (conv<FT>(v[k]) - conv<FT>(c[k])) * (conv<FT>(v[k]) - conv<FT>(c[k])) * bstarnorm[k];
		//cout << v << "@" << k << " "<< rho << " / " << c[k] << endl;
		if(rho[k] < R) {
			if(k==0)
			{
				R = rho[0];
				v_min = v;
				cout << "R: " << R << endl;
				if (R==3) {
					cout << "Returning." << endl;
					break;
				}
			}
			else {
				k--;
				//r[k-1] = max(r[k-1], r[k]); // From Paper

				for(int j = r[k + 1]; j > k; j--) {
					sigma[k][j] = sigma[k][j+1] - conv<FT>(v[j]) * mu[j][k];
				}
				r[k]     = max(r[k], r[k + 1]);
				r[k + 1] = k + 1;


				c[k] = sigma[k][k+1];  // -muProdDouble<ZT, FT>(v, mu, k, dim-1);
				v[k] = conv<FT>(round(c[k]));
				w[k] = 1;
			}
		}
		else {
			k++;
			if(k == dim) {
				cerr << "There is no solution.\n";
				break;
			}
			r[k-1] = k;
			if (k >= last_nonzero) {
				last_nonzero = k;
				v[k] = v[k] + 1;
			}
			else {
				if(conv<FT>(v[k]) > c[k]) {
					v[k] = v[k] - w[k];
				}
				else {
					v[k] = v[k] + w[k];
				}
				w[k] = w[k] + 1;
			}
		}
	} // End while

	calcVectorLengthDouble<ZT, FT>(v_min, B);
}

template<class ZT, class FT> ZT next(ZT a, FT r) {
	ZT am1 = a - 1;
	ZT ap1 = a + 1;
	ZT ap;
	FT left = abs(conv<FT>(am1) - r);
	FT right = abs(conv<FT>(ap1) - r);
	FT middle = abs(conv<FT>(a) - r);

	if (left < right)
		ap = am1;

	else
		ap = ap1;

	//if (middle < left && middle < right)
	//	ap = a;

	if(!     (    (a<r) && (r<ap)      )) {
			cerr << "Condition" << endl;
			cout << a << " " << r << " " << ap << endl;
	}

	return ap;
}

*
 * Following Parallel Enumeration of Shortest Lattice Vectors
 * Algorithm 1: Basic Enumeration Algorithm

template<class ZT, class FT> FT DagdelenEnumeration(const Mat<ZT>& B, const Mat<FT>& mu, const Vec<FT>& bstarnorm,
		Vec<ZT>& u, Vec<FT>& l, Vec<FT>& c, Vec<ZT>& v, Vec<ZT>& Delta, Vec<ZT>& delta,
		const FT& A_in, int j, int k, int dim, bool fullenum=true) {
	//Vec<ZT> u; u.SetLength(dim + 1);//(k - j + 2);

	Vec<ZT> umin; umin.SetLength(dim + 1);//(k - j + 2);
	Vec<ZT> vmin; vmin.SetLength(dim + 1);//(k - j + 2);

	Vec<FT> lmin; lmin.SetLength(dim + 1);//(k - j + 2);
	Vec<FT> cmin; cmin.SetLength(dim + 1);//(k - j + 2);

	// To decide the zigzag pattern
	Vec<ZT> Deltamin; Deltamin.SetLength(dim + 1);//(k - j + 2);
	Vec<ZT> deltamin; deltamin.SetLength(dim + 1);//(k - j + 2);

	int s, t; // s required for zigzagpattern
	s = t = j;

	umin = u;

	// Set scalars
	FT A = A_in;
	FT l_min = FT(0);

	while(t <= k) {
		l[t] = l[ t + 1 ] + ( conv<FT>(u[t]) + c[t] ) * ( conv<FT>(u[t]) + c[t] ) * bstarnorm[t];
		//cout << u << c << " " << t << endl;
		if( l[t] < A) {

			if (t > j) {
				t = t - 1;
				c[t] = muProd<ZT, FT>(u, mu, t, k);
				u[t] = v[t] = - conv<ZT>(round ( c[t] )); //

				// For zigzag
				Delta[t] = 0;
				if(conv<FT>(u[t]) > - c[t]) {
					delta[t] = - 1;
				}
				else {
					delta[t] = 1;
				}
			}

			else {
				A = l[t];
				umin = u; lmin = l; cmin = c;
				vmin = v; Deltamin = Delta; deltamin = delta;
				if(!fullenum)
					return A;
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

	u = umin;	l = lmin;	cmin = c;
	v = vmin;	Delta = Deltamin; delta = deltamin;
	return A;
}


template<class ZT, class FT> void SchnorrEnumeration(const NTL::Mat<ZT>& B, const NTL::Mat<FT>& mu, const NTL::Vec<FT>& bstarnorm, const FT& A, int j, int k, int dim, NTL::Vec<ZT>& u) {
	NTL::Vec<FT> cbar; cbar.SetLength(dim + 1);
	NTL::Vec<FT> ctilde; ctilde.SetLength(dim + 1);
	NTL::Vec<ZT> utilde; utilde.SetLength(dim + 1);
	NTL::Vec<ZT> v; v.SetLength(dim + 1);
	NTL::Vec<FT> y; y.SetLength(dim + 1);
	// To decide the zigzag pattern
	NTL::Vec<ZT> Delta; Delta.SetLength(dim + 1);
	NTL::Vec<ZT> delta; delta.SetLength(dim + 1);

	utilde[j] = u[j] = 1;
	y[j] = 0;
	v[j] = 0;
	Delta[j] = 0;
	delta[j] = 1;
	int s, t;
	s = t = j;

	cbar[j] = conv<FT>(bstarnorm[j]) * 1.001;

	for(int i=j+1; i <= dim; i++) {
		ctilde[i] = y[i] = 0;
		u[i] = utilde[i] = Delta[i] = v[i] = 0;
		delta[i] = 1;
		cbar[i] = conv<FT>(bstarnorm[j]) * 1.01;
	}

	for(int i=j; i <= k; i++) {
		cbar[i] = bstarnorm[0] * 1.001;// * min(1.0, 1.05*(dim - i) / (dim - 1));
	}

	while (t <= k) {
		ctilde[t] = ctilde[t+1] + (y[t] + conv<FT>(utilde[t])) * (y[t] + conv<FT>(utilde[t])) * conv<FT>(bstarnorm[t]);
		//cout << "u: " << u << "c: " << ctilde << "y:" << y[t] << "t: " << t << endl;
		// Actual norm is smaller than target
		if(ctilde[t] < cbar[t]) {

			// We are not at the root
			if(t > j) {
				t--; // Walk down the tree
				y[t] = muProd<ZT, FT>(utilde, mu, t, s);
				utilde[t] = v[t] = conv<ZT>(ceil(-y[t] - 0.5));

				//Delta[t] = 0;
				Delta[t] = 0;
				if(conv<RR>(utilde[t]) > -y[t]) {
					delta[t] = - 1;
				}
				else {
					delta[t] = 1;
				}
			}

			// Found a new short solution
			else {

				u = utilde;
				for(int i=j; i <= k; i++) {
					cbar[i] = ctilde[j];// * min(1.0, 1.05*(dim - i) / (dim - 1)); // set new bound
				}
			}
		}

		// Size already too large
		else {
			t++; // go up the tree
			s = max(s,t);

			if(t<s) Delta[t] = -Delta[t];
			if(Delta[t] * delta[t] >= 0) Delta[t] += delta[t];

			utilde[t] = v[t] + Delta[t];

			s = max(s,t); //

			if(t < s) {
			    if(Delta[t] == 0) Delta[t]++;
			    else if(Delta[t] < 0) Delta[t] = abs(Delta[t]) + 1;
			    else if(Delta[t] > 0) Delta[t] = - Delta[t];
			}
			else {
			    Delta[t]++;
			}
			utilde[t] = v[t] + Delta[t];
		}
	}
	cout << "A: " << cbar[0] << endl;
}*/




int main(int argc, char** argv)
{
	//searchShortestSeedCandidates(100, 1, 4000);
	/*uint32_t val=0;
	val=val | 0x0001;
	val=val | 0x0002;
	val=val | 0x0004;

	cout << val << endl;
	return 0;*/

#pragma omp parallel
#pragma omp single nowait
	cout << "Working threads: " << omp_get_num_threads() << endl;

#ifdef USE_MB_TYPES
	cout << "Using MB types." << endl;
#else
	cout << "Using NTL types." << endl;
#endif

	if(readConfig("") < 0)
		readConfig("/home/mburger/Dokumente/p3enum/penum/");
		

	std::string filepath=""; // Basisfile that is read
	std::string vecpath=""; // Vectorfile that is read

    int mode = 1;     	// 1 = Do basis reduction and search
                    	// 2 = Do LLL/BKZ and write base to file

	int candidate_height=2; // How long are the candidate vectors searched in serial


	// Parse input arguments
    for(int i = 1; i < argc; i++) {
    	std::string input = std::string(argv[i]);

        if(input == "--basisfile") {
            if(argc <= i) {
                std::cerr << "Missing basisfile in argument. Terminating." << std::endl;
                exit(-1);
            }

            filepath = std::string(argv[i+1]);
            i++;
        }

        if(input == "--prunefile") {
            if(argc <= i) {
                std::cerr << "Missing prunefile in argument. Terminating." << std::endl;
                exit(-1);
            }

            Configurator::getInstance().ext_pruning_function_file = std::string(argv[i+1]);
            i++;
        }

        if(input == "--vecfile") {
            if(argc <= i) {
                std::cerr << "Missing vectorfile in argument. Terminating." << std::endl;
                exit(-1);
            }

            cout << "Only executing vector length mode." << endl;
            vecpath = std::string(argv[i+1]);
            i++;
        }

        if(input == "--beta") {
            if(argc <= i) {
                std::cerr << "Missing beta in argument. Terminating." << std::endl;
                exit(-1);
            }

            int beta = atoi(argv[i+1]);
            Configurator::getInstance().glob_beta = beta;
            i++;
        }

        if(input == "--sheight") {
            if(argc <= i+1) {
                std::cerr << "Missing sheight in argument. Terminating." << std::endl;
                exit(-1);
            }

            candidate_height = atoi(argv[i+1]);
            Configurator::getInstance().forced_height = candidate_height;
            i++;
        }
        
        if(input == "--prebeta") {
             if(argc <= i+1) {
                std::cerr << "Missing prebeta in argument. Terminating." << std::endl;
                exit(-1);
            }   

            int prebeta = atoi(argv[i+1]);
            Configurator::getInstance().prebeta = prebeta;
            cout << "Setting prebeta parameter to " << prebeta << "." << endl;
            i++;  
        }       
        
        if(input == "--pruneparam") {
             if(argc <= i+1) {
                std::cerr << "Missing pruneparam in argument. Terminating." << std::endl;
                exit(-1);
            }   

            double prune_param = atof(argv[i+1]);
            Configurator::getInstance().prune_param = prune_param;
            cout << "Setting pruning parameter to " << prune_param << "." << endl;
            i++;  
        }

        if(input == "--enumpruning") {
            if(argc <= i+1) {
                std::cerr << "Missing enumpruning in argument. Terminating." << std::endl;
                exit(-1);
            }
            std::string min_prd(argv[i+1]);
            Configurator::getInstance().enum_prune = to_prune(min_prd);
            cout << "Activating " << min_prd << " pruning." << endl;
            i++;
        }

        if(input == "--nogaussestimate") {
        	Configurator::getInstance().do_gauss_estimate = false;
        }

        if(input == "--Enumthreads") {
            if(argc <= i+1) {
                std::cerr << "Missing Enumthreads in argument. Terminating." << std::endl;
                exit(-1);
            }

            Configurator::getInstance().Enumthreads = atoi(argv[i+1]);
            i++;
        }        
        
        if(input == "--candsize") {
            if(argc <= i+1) {
                std::cerr << "Missing candsize in argument. Terminating." << std::endl;
                exit(-1);
            }

            Configurator::getInstance().cand_queue_size = atoi(argv[i+1]);
            i++;
        }

        if(input == "--parthres") {
            if(argc <= i+1) {
                std::cerr << "Missing parthres in argument. Terminating." << std::endl;
                exit(-1);
            }

            Configurator::getInstance().par_threshold = atoi(argv[i+1]);
            i++;
        }

        if(input == "--delta") {
            if(argc <= i) {
                std::cerr << "Missing delta in argument. Terminating." << std::endl;
                exit(-1);
            }

            double delta = atof(argv[i+1]);
            Configurator::getInstance().glob_delta = delta;
            i++;
        }

        if(input == "--amax") {
            if(argc <= i) {
                std::cerr << "Missing amax in argument. Terminating." << std::endl;
                exit(-1);
            }

            double amax = atof(argv[i+1]);
            Configurator::getInstance().Amax = amax;
            i++;
        }


        if(input == "--aend") {
            if(argc <= i) {
                std::cerr << "Missing aend in argument. Terminating." << std::endl;
                exit(-1);
            }

            double aend = atof(argv[i+1]);
            Configurator::getInstance().Aend = aend;
            i++;
        }

        if(input == "--mbenum") {
            Configurator::getInstance().use_ntl_enum = false;
        }

        if(input == "--itenum") {
        	Configurator::getInstance().iterative_enumeration = true;
        }

        if(input == "--onlyLLLexport") {
            if(argc <= i) {
                std::cerr << "Missing output file argument for LLL. Terminating." << std::endl;
                exit(-1);
            }

            cout << "Base will be exported to: " << std::string(argv[i+1]) << endl;
            Configurator::getInstance().reducedfilepath = std::string(argv[i+1]);
            mode = 2;
            i++;
        }

        if(input == "--noLLL") {
        	Configurator::getInstance().dolll = false;
        }
    }

    std::cout << "Delta: " << Configurator::getInstance().glob_delta << std::endl;
    std::cout << "Beta: " << Configurator::getInstance().glob_beta << std::endl;
    std::cout << "Mode: " << mode << std::endl;
    std::cout << "Candidate Height: " << candidate_height << std::endl;
    std::cout << "Target length: " << Configurator::getInstance().Amax << std::endl;
    std::cout << "No. of threads in ENUM: " << Configurator::getInstance().Enumthreads << std::endl;

    // Read lattice basis from file
    if(filepath == "") {
    	std::cerr << "No valid filepath submitted. Exiting." << std::endl;
    	exit(-1);
    }

    // Execute if the parallel randomizing approach ist chosen


    // Standard procedure
    NTL::mat_ZZ B;
    std::ifstream basestream;
    std::cout << "Trying to read: " << filepath << std::endl;
    basestream.open(filepath.c_str(), ios_base::in);
    try {
    	basestream >> B;
    }
    catch (std::exception & e) {
    	cout << "Error while reading basis. Exiting" << std::endl;
    	exit(-1);
    }
    basestream.close();

    std::cout << "Read a basis of dimension: "
    		<< B.NumRows() << " / " << B.NumCols()
			<< " with |b_0| = " << lengthOfVector<ZZ, RR> (B[0])
			<< std::endl;

    if(vecpath!="") {
    	NTL::vec_ZZ resvec;
    	std::ifstream vecstream;
    	vecstream.open(vecpath.c_str(), ios_base::in);
    	cout << "Reading vec: " << vecpath << endl;
        try {
        	vecstream >> resvec;
        	cout << resvec << endl;;
        }
        catch (std::exception & e) {
        	cout << "Error while solution vector. Exiting" << std::endl;
        	exit(-1);
        }
        cout << "Length of input vector: " << lengthOfVector<ZZ,RR>(resvec) << endl;;
        exit(0);
    }

    //LLL_XD(B, Configurator::getInstance().glob_delta);
    cout << "Starting procedure" << endl;

	int dim = B.NumCols();

	pEnumeratorDouble pener = pEnumeratorDouble(dim);

	vec_ZZ vec; vec.SetLength(B.NumRows() + 3);

	pener.solveSVP(B,vec);


   return 0;
}
