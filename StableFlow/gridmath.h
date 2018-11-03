#pragma once
#include "shared.h"

Float Grid_Norm(const Grid &G);

enum PreConditioner { DIAG };

template<typename T, PreConditioner prec>
class MFPCG {
public:
	// See: Fluid simulation for computer graphics, Robert Bridson, Chapter 4.3
	int n, m;
	int size;
	T *Adiag;
	T *Aplusi;
	T *Aplusj;
	T *b;//RHS constant
	T *p;//guess of answer
	T *s;//search vector
	T *z;//auxillary vector
	T *r;//residual
	MFPCG();
	~MFPCG();
	int idx(int i, int j);
	T Dot_Product(T * A, T * B);
	T Norm_Inf(T *x);
	void Apply_A_To(T * x, T * ax);
	void Apply_Prec_To(T * x, T * px);
	void Vector_Comb_To(T * a, T * b, T c, T * x);//x=a+b*c;
	void solve(Float tolerance);
};

template<typename T>
void Load_Array_From_Grid(T *A, const Grid &B) {
	int sz = B.size();
	for (int i = 0; i < sz; i++) {
		A[i] = B(i);
	}
}

template<typename T, PreConditioner prec>
int MFPCG<T, prec>::idx(int i, int j) {
	return i * m + j;
}

template<typename T, PreConditioner prec>
T MFPCG<T, prec>::Dot_Product(T *A, T *B) {
	T ret = 0;
	for (int i = 0; i < size; i++) {
		ret += A[i] * B[i];
	}
	return ret;
}

template<typename T, PreConditioner prec>
T MFPCG<T, prec>::Norm_Inf(T * x) {
	T ret = 0;
	for (int i = 0; i < size; i++) {
		ret = max(ret, fabs(x[i]));
	}
	return ret;
}

template<typename T, PreConditioner prec>
void MFPCG<T, prec>::Apply_A_To(T* x, T* ax) {
	// A*x = ax
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			int t = idx(i, j);
			ax[t] = 0;
			//A[(i,j),(i-1,j)]
			if (i > 0) ax[t] += Aplusi[idx(i - 1, j)] * x[idx(i - 1, j)];
			//A[(i,j),(i,j-1)]
			if (j > 0) ax[t] += Aplusj[idx(i, j - 1)] * x[idx(i, j - 1)];
			//A[(i,j),(i,j)]
			ax[t] += Adiag[idx(i, j)] * x[idx(i, j)];
			//A[(i,j),(i+1,j)]
			if (i + 1 < n) ax[t] += Aplusi[idx(i, j)] * x[idx(i + 1, j)];
			//A[(i,j),(i,j+1)]
			if (j + 1 < m) ax[t] += Aplusj[idx(i, j)] * x[idx(i, j + 1)];
		}
	}
}

template<typename T, PreConditioner prec>
void MFPCG<T, prec>::Apply_Prec_To(T * x, T * px) {
	if (prec == DIAG) {
		for (int i = 0; i < size; i++) {
			px[i] = x[i] * Adiag[i];
		}
	}
}

template<typename T, PreConditioner prec>
void MFPCG<T, prec>::Vector_Comb_To(T * a, T * b, T c, T * x) {
	for (int i = 0; i < size; i++) {
		x[i] = a[i] + b[i] * c;
	}
}

template<typename T, PreConditioner prec>
MFPCG<T, prec>::MFPCG() {
	n = MESHH - 2;
	m = MESHW - 2;
	size = n * m;
	int ary_size = sizeof(T)*size;
	Adiag = (T*)malloc(ary_size);
	Aplusi = (T*)malloc(ary_size);
	Aplusj = (T*)malloc(ary_size);
	b = (T*)malloc(ary_size);
	p = (T*)malloc(ary_size);
	s = (T*)malloc(ary_size);
	z = (T*)malloc(ary_size);
	r = (T*)malloc(ary_size);
}

template<typename T, PreConditioner prec>
MFPCG<T, prec>::~MFPCG() {
	free(Adiag);
	free(Aplusi);
	free(Aplusj);
	free(b);
	free(p);
	free(s);
	free(z);
	free(r);
}

template<typename T, PreConditioner prec>
void MFPCG<T, prec>::solve(Float tolerance) {
	SpMat MA(size, size);
	VectorXd Mb(size);
	for (int i = 0; i < size; i++) Mb(i) = b[i];
	vector<Tri> v;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			int eid = idx(i, j);
			v.push_back(Tri(eid, eid, Adiag[eid]));
			if (i + 1 < n) {
				int eid1 = idx(i + 1, j);
				v.push_back(Tri(eid, eid1, Aplusi[eid]));
			}
			if (j + 1 < m) {
				int eid1 = idx(i, j + 1);
				v.push_back(Tri(eid, eid1, Aplusj[eid]));
			}
		}
	}
	MA.setFromTriplets(v.begin(), v.end());
	ConjugateGradient<SparseMatrix<Float>, Upper > pcg;
	pcg.compute(MA);
	VectorXd X = pcg.solve(Mb);


	for (int i = 0; i < size; i++) p[i] = 0;//initial guess p=0
	for (int i = 0; i < size; i++) r[i] = b[i];//residual vector r=b
	Apply_Prec_To(r, z);//z=applyPreconditioner(r)
	for (int i = 0; i < size; i++) s[i] = z[i];//s=z
	T sigma = Dot_Product(z, r);
	for (int iter = 0; iter < MAX_PCG_STEP; iter++) {
		//for (int i = 0; i < 50; i++) { cout << p[i] << " "; }cout << endl;
		Apply_A_To(s, z);//z=applyA(s)
		T alpha = sigma / Dot_Product(z, s);
		Vector_Comb_To(p, s, alpha, p);//p=p+alpha*s
		//for (int i = 0; i < 50; i++) { cout << s[i] << " "; }cout << endl;
		Vector_Comb_To(r, z, -alpha, r);
		if (Norm_Inf(r) <= tolerance) {
			break;
		}
		Apply_Prec_To(r, z);
		T sigma1 = Dot_Product(z, r);
		T beta = sigma1 / sigma;
		Vector_Comb_To(z, s, beta, s);
		sigma = sigma1;
	}
	cout << "PCG run done, r=" << Norm_Inf(r) << endl;
	// The calculated answer is p

	VectorXd xp(size);
	for (int i = 0; i < size; i++) xp(i) = p[i];
	cout << (xp - X).cwiseAbs().maxCoeff() << endl;
}

class ConstMask {
public:
	Grid dlt;
	Eigen::ArrayXXi msk;
	ConstMask() {}
	void resize(const Grid &A);
	void Set_Box(Float x0, Float x1, Float y0, Float y1, Float d);
	void Set_Ellipse(Float x0, Float y0, Float rh, Float rw, Float d);
	void Set_Real_Circle(Float x0, Float y0, Float rh, Float d);
	//ConstBlock(const Grid &A, Float x0, Float x1, Float y0, Float y1, Float _d);
};

void Apply_ConstMask(Grid &A, const ConstMask &B);
//void Add_Block(Grid &A, Float x0, Float x1, Float y0, Float y1, Float d);

void Truncate_Position(const Grid &A, Float &x, Float &y);
Float Interpolate(const Grid &A, Float x, Float y);

//For now we use a periodical boundary condition
Grid Grid_D(const Grid &A, AXIS ax, DIFFTYPE typ);
Grid Grid_S(const Grid &U, const Grid &V, const Grid &phi);
//Float Local_Diff(Matrix A, int x, int y, int order, AXIS ax, DIFFTYPE typ);

