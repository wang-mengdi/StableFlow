#pragma once
#include "shared.h"
#include "gridmath.h"
//#include "fiber.h"

typedef Eigen::Triplet<Float> Tri;
typedef Eigen::SparseMatrix<Float> SpMat;

class Dye {
public:
	Grid dens;
	Grid src;
	Float rgb[3];
	void Init(int _r, int _g, int _b);
};

class Solver {
public:
	//Currently we use traditional, not staggered mesh.
	Grid U, V;//x/y component of velocity
	Grid U1, V1;
	Grid P,div;
	Grid FU, FV;
	vector<Dye> colors;
	int step_time=0;
	void Init(void);
	void Draw_Colors(int screenid);
	void Draw_Velocity_Field(const Grid & U, const Grid & V, int screenid);
	void Draw(void);
	void Step();
	void Diffuse(AXIS ax, Grid & x, const Grid & b, Float diff);
	void Velocity_Step(void);
	void Density_Step(Grid & D, const Grid & S);
	void Step_Fluid(void);
};