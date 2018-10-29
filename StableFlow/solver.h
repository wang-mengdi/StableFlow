#pragma once
#include "shared.h"
#include "gridmath.h"
//#include "fiber.h"

typedef Eigen::Triplet<Float> Tri;
typedef Eigen::SparseMatrix<Float> SpMat;


class Dye {
public:
	Grid dens;
	ConstMask src;
	Float rgb[3];
	Dye() {}
	Dye(int _r, int _g, int _b);
};

class Solver {
public:
	//Currently we use traditional, not staggered mesh.
	MFPCG<Float, DIAG> pcg;
	Grid U, V;//x/y component of velocity
	Grid U1, V1;
	Grid P,div;
	vector<Dye> colors;
	ConstMask CU;
	ConstMask CV;
	int step_time=0;
	void Init(void);
	void Draw_Colors(int screenid);
	void Draw_Velocity_Field(const Grid & U, const Grid & V, int screenid);
	void Draw(void);
	void Step();
	void Diffuse(AXIS ax, Grid & x, const Grid & b, Float diff);
	void Load_Pressure_System(const Grid & div);
	void Fill_Pressure_Solution(Grid & P);
	void Get_Pressure(Grid & P, const Grid & U, const Grid & V, Grid & div);
	void Project(Grid & U, Grid & V, Grid & P, Grid & div);
	void Velocity_Step(void);
	void Density_Step(Dye & D);
	void Step_Fluid(void);
};