#pragma once
#include "shared.h"

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