#include "p3BasisReduction.hpp"

using namespace p3;
using namespace std;

/**
 * Calculates SS(B) as defined in [1], Page 2490, after Theorem 1
 */
template<class T, class FT>
FT p3BasisReduction<T, FT>::calc_SS_B (const GSInfo<T, FT>& gs_info) {
	FT res = FT(0);
	for (unsigned int i = 0; i < gs_info._bstarabs.size(); i++) {
		res += gs_info._bstarabs[i];
	}

	return res;
}

/**
 * Calculates D_l^(k) as defined in [1], Page 2493, Theorem 3
 */
template<class T, class FT>
FT p3BasisReduction<T, FT>::calc_D_lk(p3Matrix<FT>& fmat, GSInfo<T,FT>& gsinfo, int l, int k) {
	FT res = FT(0);
	for(int j=l; j<=k; j++) {
		res += (gsinfo._mu[k][j] * gsinfo._mu[k][j]) * (gsinfo._bstarabs[j]);
	}
	return res;
}

template<class T, class FT>
	void p3BasisReduction<T, FT>::afterDeepInsertUpdate (const int i, const int k, const int n, GSInfo<T, FT>& gsinfo) {

	FT tmp1;

	p3Vector<FT> D, Dinv, P, epsilons, S;
	D.SetLength(n); Dinv.SetLength(n);
	P.SetLength(n);
	epsilons.SetLength(n);
	S.SetLength(n);

	P[k] = gsinfo._bstarabs[k]; // P_k = B_k
	D[k] = gsinfo._bstarabs[k]; // D_k = B_k

	// 2:4
	for (int j = k - 1; j >= i; --j) {
		P[j] = gsinfo._mu[k][j] * gsinfo._bstarabs[j];
		D[j] = D[j + 1] + gsinfo._mu[k][j] * P[j];
	}


	// 5 Set S_i = S_i+1 = ... = S_n = 0
	for (int c = i; c < n; c++) {
		S[c] = 0;
	}

	// 30:32
	Dinv[k] = 1.0 / D[k]; // da mu[k][k] == 1

	for (int j = k; j >= i + 1; --j) {
		Dinv[j-1] = 1.0 / D[j-1];
		gsinfo._bstarabs[j] = (D[j] * gsinfo._bstarabs[j-1]) * Dinv[j-1];
	}
	gsinfo._bstarabs[i] = D[i];

	// 8
	tmp1 = gsinfo._mu[k][k-1] * Dinv[k];
	// 9:11
	for (int l = n - 1; l >= k+1; --l) {
		S[l] = P[k] * gsinfo._mu[l][k];
		gsinfo._mu[l][k] = gsinfo._mu[l][k-1] - tmp1 * S[l];
	}

	// 7
	for (int j = k - 1; j >= i+1; --j) {
		// 8
		tmp1 = gsinfo._mu[k][j-1] * Dinv[j];

		// 9:11
		for (int l = k + 1; l < n; ++l) {
			S[l] += P[j] * gsinfo._mu[l][j];
			gsinfo._mu[l][j] = gsinfo._mu[l][j-1] - tmp1 * S[l];
		}

		// 12:14
		for (int l = k; l >= j + 2; --l) {
			S[l] += P[l] * gsinfo._mu[l-1][j];
			gsinfo._mu[l][j] = gsinfo._mu[l-1][j-1] - tmp1 * S[l];
		}

		// 15:17
		S[j+1] = P[j];
		gsinfo._mu[j+1][j] = gsinfo._mu[j][j-1] - tmp1 * S[j+1];
	}

	// 21
	for (int l = k + 1; l < n; ++l) {
		gsinfo._mu[l][i] = (S[l] + P[i] * gsinfo._mu[l][i]) * Dinv[i];
	}

	// 22
	for (int l = k; l >= i+2; --l) {
		gsinfo._mu[l][i] = (S[l] + P[i] * gsinfo._mu[l-1][i]) * Dinv[i];
	}

	// 23
	gsinfo._mu[i+1][i] = P[i] * Dinv[i];

	// Copy von Zeile 26
	for(int j = 1; j <= i-1; ++j) {
		epsilons[j] = gsinfo._mu[k][j];
	}

	// 25:29
	for(int j = k; j != i; --j) {
		for(int l = 1; l < i; ++l) {
			gsinfo._mu[j][l] = gsinfo._mu[j-1][l];
		}
	}

	// Copy von Zeile 28
	for(int j = 1; j < i; ++j) {
		gsinfo._mu[i][j] = epsilons[j];
	}

}

template<class T, class FT>
	void p3BasisReduction<T, FT>::reduce_RET(const int k, const int l, p3Matrix<T>& mat, p3Matrix<FT>& fmat, p3Matrix<T>& H,  GSInfo<T, FT>& gs_info) {

	if( abs(gs_info._mu[k][l]) <= 0.50001)
		return;

	FT q_f = round( gs_info._mu[k][l] );

	fmat[k] = fmat[k] - q_f * fmat[l];
	H[k] = H[k] - q_f * H[l];
	gs_info._mu[k][l] = gs_info._mu[k][l] - q_f;

	for (int i = 0; i < l; i++) {
		gs_info._mu[k][i] = gs_info._mu[k][i] - q_f*gs_info._mu[l][i];
	}

	return;
}

template<class T, class FT>
	void p3BasisReduction<T, FT>::swap_SWAP(const int k, const int k_max, p3Matrix<T>& mat, p3Matrix<FT>& fmat, p3Matrix<T>& H,  GSInfo<T, FT>& gs_info) {
	assert( k > 0);

	// Exchange b_k and b_k-1
	p3Vector<FT> tvec_f = fmat[k];
	fmat[k] = fmat[k-1];
	fmat[k-1] = tvec_f;

	// Exchange H[k] and H[k-1]
	p3Vector<T> tvec_i = H[k];
	H[k] = H[k-1];
	H[k-1] = tvec_i;

	if(k>1)
	{
		for(int j=0; j <= k-2; j++) {
			FT tmu_f = gs_info._mu[k][j];
			gs_info._mu[k][j] = gs_info._mu[k-1][j];
			gs_info._mu[k-1][j] = tmu_f;
		}
	}

	FT mu = gs_info._mu[k][k-1];
	FT B = gs_info._bstarabs[k] + mu * mu * gs_info._bstarabs[k-1];
	gs_info._mu[k][k-1] = mu * gs_info._bstarabs[k-1] / B;

	p3Vector<FT> b = gs_info._bstar[k-1];
	gs_info._bstar[k-1] = gs_info._bstar[k] + mu*b;
	gs_info._bstar[k] = - gs_info._mu[k][k-1] * gs_info._bstar[k] + (gs_info._bstarabs[k] / B) * b;
	gs_info._bstarabs[k] = gs_info._bstarabs[k-1] * gs_info._bstarabs[k] / B;
	gs_info._bstarabs[k-1] = B;

	for (int i = k+1; i <= k_max; i++) {
		FT t = gs_info._mu[i][k];
		gs_info._mu[i][k] = gs_info._mu[i][k-1] - mu * t;
		gs_info._mu[i][k-1] = t + gs_info._mu[k][k-1] * gs_info._mu[i][k];
	}

	return;
}

template<class T, class FT>
	void p3BasisReduction<T, FT>::sizeReduce(const int k, p3Matrix<T>& mat, p3Matrix<FT>& fmat, GSInfo<T, FT>& gs_info) {

	for(int j =  k-1; j >= 0; j--) {

		if(fabs(gs_info._mu[k][j]) > 0.5) {
			FT q_j = round( gs_info._mu[k][j] );

			// RR to ZZ
			//std::string strval = convert2Str<FT>(q_j);
			//T q_j_i; convertFromStr<T>(q_j_i, strval);
			//mat[k] = mat[k] - q_j_i * mat[j];

			fmat[k] = fmat[k] - q_j * fmat[j];
			//_fmat = conv_to_ft_mat <T, FT> (_mat);

			// Update only relevant mu
		    //for (int i = 0; i < j; i++)
		    //	gs_info._mu[k][i] = gs_info._mu[k][i] - q_j * gs_info._mu[j][i];
			GramSchmidt(mat, fmat, gs_info);
		}
		// Full update
		GramSchmidt(mat, fmat, gs_info);
	}
}

template<class T, class FT>
	void p3BasisReduction<T, FT>::GramSchmidt (p3Matrix<T>& mat, p3Matrix<FT>& fmat, GSInfo<T, FT>& gs_info) {

	gs_info._bstar[0] = fmat[0];
	for(int i=1; i < (int)fmat.NumRows(); i++) {
		p3Vector<FT> v = fmat[i];

		for(int j=i-1; j>= 0; j--) {
			gs_info._mu[i][j] = (fmat[i] * gs_info._bstar[j]) / ( gs_info._bstar[j] * gs_info._bstar[j]);
			v = v - (gs_info._mu[i][j] * gs_info._bstar[j]);
		}
		gs_info._bstar[i] = v;
	}

	for(int i=0; i < (int)fmat.NumRows(); i++) {
		gs_info._mu[i][i] = 1.0;
	}

	// Fill ||b*||^2
	assert(gs_info._bstarabs.length() == gs_info._bstar.NumRows());

	for(unsigned int c=0; c < gs_info._bstar.NumCols(); c++) {
		gs_info._bstarabs[c] = gs_info._bstar[c] * gs_info._bstar[c];
	}

	return;
}

template<class T, class FT>
void p3BasisReduction<T, FT>::hello(T i)
{
	cout << "Hello" << T(i)  << endl;
}

template<class T, class FT>
int p3BasisReduction<T, FT>::sigma_i_k(int i, int k, p3Matrix<T>& mat, p3Matrix<FT>& fmat) {

	p3Vector<FT> row_back = fmat[i];
	p3Vector<FT> row_back2 = fmat[i];
	p3Vector<FT> row_k = fmat[k];

	fmat[i] = row_k;

	for(int pos=i; pos < k; pos++) {
		row_back2 = fmat[pos+1];
		fmat[pos+1] = row_back;
		row_back = row_back2;
	}


	/*p3Vector<T> row_backi = _mat[i];
	p3Vector<T> row_back2i = _mat[i];

	p3Vector<T> row_ki = _mat[k];

	_mat[i] = std::move(row_ki);

	for(int pos=i; pos < k; pos++) {
		row_back2i = _mat[pos+1];
		_mat[pos+1] = row_backi;
		row_backi = row_back2i;
	}*/

	return 0;
}

/**
 *
 */
// (p3Matrix<FT>& fmat, GSInfo<T, FT>& gs_info, std::vector<FT>& D_j_k, int i, int k);
template<class T, class FT>
FT p3BasisReduction<T, FT>::calc_S_ik(p3Matrix<FT>& fmat, GSInfo<T, FT>& gsinfo, std::vector<FT>& D_j_k, int i, int k) {
	FT ret = FT(0);

	for(int j=i; j<=k-1; j++) {
		ret += gsinfo._mu[k][j] * gsinfo._mu[k][j] * gsinfo._bstarabs[j] * (gsinfo._bstarabs[j]/D_j_k[j] - 1);
	}

	return ret;
}


/**
 *
 */
template<class T, class FT>
int p3BasisReduction<T, FT>::argmax_S_lk(p3Matrix<FT>& fmat, GSInfo<T, FT>& gs_info, int k, std::vector<FT>& S_lk, int& i) {
	i = -1;

	// Calculate all S_lk with 1 <= l <= k-1

	// 1. Intermediate vector with D_lk values
	std::vector<FT> D_lk;
	for(int l=0; l <= k; l++) {
		D_lk.push_back(calc_D_lk(fmat, gs_info, l, k));
	}

	// 2. Fill the vector with S_lk-values
	for(int l=0; l <= k-1; l++) {
		S_lk.push_back(calc_S_ik(fmat, gs_info, D_lk, l, k));
	}

	// get the index of the maximum entry
	auto max_it = std::max_element(S_lk.begin(), S_lk.end());
	i = std::distance(S_lk.begin(),max_it);

	assert (i >= 0);
	assert (S_lk.size() > 0);

	return 0;
}

template<class T, class FT>
void p3BasisReduction<T, FT>::swapBaseVectors (p3Matrix<FT>& fmat, int idx1, int idx2) {
	p3Vector<FT> vec1 = fmat[idx1];
	p3Vector<FT> vec2 = fmat[idx2];
	fmat[idx1] = vec2;
	fmat[idx2] = vec1;

	//Vec<T> vec1i = _mat[idx1];
	//Vec<T> vec2i = _mat[idx2];
	//_mat[idx1] = vec2i;
	//_mat[idx2] = vec1i;

	return;
}


template<class T, class FT>
	void p3BasisReduction<T, FT>::S2LLL (p3Matrix<T>& mat, double eta) {

	assert(eta > 0.6 && eta <= 1);

	p3Matrix<FT> fmat = conv_to_ft_mat <T, FT> (mat);
	GSInfo<T, FT> gsinfo = GSInfo<T, FT>(mat.NumRows(), mat.NumCols());

	GramSchmidt(mat, fmat, gsinfo);

	FT startBB = calc_SS_B(gsinfo);
	cout << startBB << endl;
	int redcnt = 0;

	int k = 1;
	while (k < (int) mat.NumRows() ) {
		sizeReduce(k, mat, fmat, gsinfo);
		redcnt++;

		// Compute i <-argmax
		int i = -1;
		std::vector<FT> S_jk;
		argmax_S_lk(fmat, gsinfo, k, S_jk, i);

		if (S_jk[i] <= (1-eta) * calc_SS_B(gsinfo)) {
			k++;
		}

		else {
			sigma_i_k(i,k, mat, fmat);
			afterDeepInsertUpdate(i, k, mat.NumRows(), gsinfo);
			k = std::max(i, 1);
		}
	}

	// Copy results to the integer-representation
	mat = conv_to_int_mat<FT, T> (fmat);

	//std::cout << fmat << endl;
	cout << gsinfo._bstarabs;
	cout << "Det End:" << gsinfo.det() << endl;
	FT endBB = calc_SS_B(gsinfo);
	std::cout << "Improved to: " << endBB << " / " << endBB/startBB << std::endl;
	std::cout << redcnt << " reductions." << std::endl;
}

template<class T, class FT>
	void p3BasisReduction<T, FT>::LLL_FP (p3Matrix<T>& mat, double delta_) {
	FT delta = FT(delta_);

	GSInfo<T, FT> gsinfo = GSInfo<T, FT>(mat.NumRows(), mat.NumRows());

	p3Matrix<FT> fmat = conv_to_ft_mat <T, FT> (mat);
	GramSchmidt(mat, fmat, gsinfo);
	int k=1;
	int redcnt = 0;
	FT startBB = calc_SS_B(gsinfo);
	cout << "Start SS " << startBB << endl;
	while (k < (int)mat.getNumRows()) {
		sizeReduce(k, mat, fmat, gsinfo);
		redcnt++;

		if (gsinfo._bstarabs[k] >= (delta - gsinfo._mu[k][k-1] * gsinfo._mu[k][k-1]) * gsinfo._bstarabs[k-1] ) {
			k = k + 1;
		}

		else {
			swapBaseVectors(fmat, k, k-1);
			calc_SS_B(gsinfo);
			k = std::max(1,k-1);
		}
	}

	cout << fmat << endl;
	FT endBB = calc_SS_B(gsinfo);
	std::cout << "Improved to: " << endBB << " / " << endBB/startBB << std::endl;
	std::cout << redcnt << " reductions." << std::endl;
}

template<class T, class FT>
	void p3BasisReduction<T, FT>::LLL_FP_COHEN (p3Matrix<T>& mat, FT delta_) {
	p3Matrix<FT> fmat = conv_to_ft_mat <T, FT> (mat);
	GSInfo<T, FT> gsinfo = GSInfo<T, FT>(mat.NumRows(), mat.NumRows());

	GramSchmidt(mat, fmat, gsinfo);
	cout << gsinfo._bstar << endl;
	cout << "Det Start:" << gsinfo.det() << endl;

	// Initialize
	int k = 1;
	int k_max = 0;
	gsinfo._bstar[0] = fmat[0];
	gsinfo._bstarabs[0] = fmat[0] * fmat[0];
	p3Matrix<T>  H = p3Matrix<T>(fmat.NumRows(), fmat.NumCols(), true);
	int n = fmat.NumRows();

	while (k < n) {
		// Optional Step 2: Incremental Gram-Schmidt
		if (k > k_max) {
			k_max = k;
			gsinfo._bstar[k] = fmat[k];

			for (int j=0; j < k; j++) {
				gsinfo._mu[k][j] = fmat[k] * gsinfo._bstar[j] / gsinfo._bstarabs[j];
				gsinfo._bstar[k] = gsinfo._bstar[k] - gsinfo._mu[k][j] * gsinfo._bstar[j];
			}
			gsinfo._bstarabs[k] = gsinfo._bstar[k] * gsinfo._bstar[k];

			if(gsinfo._bstarabs[k] == 0) {
				cerr << "b_i does not form a basis" << endl;
				return;
			}
		}

		// Step 3: Test LLL-condition
		while (true) {
			reduce_RET(k, k-1, mat, fmat, H, gsinfo);

			if(gsinfo._bstarabs[k] < (delta_ - gsinfo._mu[k][k-1] * gsinfo._mu[k][k-1]) *  gsinfo._bstarabs[k-1])
			{
				// SWAP
				swap_SWAP(k, k_max, mat, fmat, H,  gsinfo);
				k = max(1, k-1);
			}

			else {
				for (int l=k-2; l >= 0; l--) {
					reduce_RET(k, l, mat, fmat, H, gsinfo);
				}
				k = k+1;
				break;
			}
		}
	}

	cout << "Det End:" << gsinfo.det() << endl;
	mat = conv_to_int_mat <FT, T> (fmat);
	return;

}

