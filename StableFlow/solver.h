#pragma once
#include "shared.h"
#include "gridmath.h"
//#include "fiber.h"

typedef Eigen::Triplet<Float> Tri;
typedef Eigen::SparseMatrix<Float> SpMat;



class Dye {
public:
	Grid dens;
	vector<ConstBlock> srcs;
	Float rgb[3];
	Dye() {}
	Dye(int _r, int _g, int _b);
	void Add_Src(Float x0, Float x1, Float y0, Float y1, Float _d);
};

class Solver {
public:
	//Currently we use traditional, not staggered mesh.
	Grid U, V;//x/y component of velocity
	Grid U1, V1;
	Grid P,div;
	vector<Dye> colors;
	vector<ConstBlock> CUs;
	vector<ConstBlock> CVs;
	int step_time=0;
	void Init(void);
	void Add_CU(Float x0, Float x1, Float y0, Float y1, Float _d);
	void Add_CV(Float x0, Float x1, Float y0, Float y1, Float _d);
	void Draw_Colors(int screenid);
	void Draw_Velocity_Field(const Grid & U, const Grid & V, int screenid);
	void Draw(void);
	void Step();
	void Diffuse(AXIS ax, Grid & x, const Grid & b, Float diff);
	void Velocity_Step(void);
	void Density_Step(Dye & D);
	void Step_Fluid(void);
};