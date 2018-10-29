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
	void PCG_run(T tolerance);
};

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