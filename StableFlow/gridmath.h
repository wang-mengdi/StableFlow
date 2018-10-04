#pragma once
#include "shared.h"

void Add_Block(Grid &A, Float x0, Float x1, Float y0, Float y1, Float d);

void Truncate_Index(const Grid &A, int &i, int &j);
Float Interpolate(const Grid &A, Float x, Float y);

//For now we use a periodical boundary condition
Grid Grid_D(const Grid &A, AXIS ax, DIFFTYPE typ);
Grid Grid_S(const Grid &U, const Grid &V, const Grid &phi);
//Float Local_Diff(Matrix A, int x, int y, int order, AXIS ax, DIFFTYPE typ);