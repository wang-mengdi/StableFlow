#pragma once
#include "shared.h"
#include "gridmath.h"
//#include "fiber.h"

typedef Eigen::Triplet<Float> Tri;
typedef Eigen::SparseMatrix<Float> SpMat;

class Solver {
public:
	//Currently we use traditional, not staggered mesh.
	Grid U, V;//x/y component of velocity
	Grid U1, V1;
	Grid P;
	Grid FU, FV;
	int step_time=0;
	void Init(void);
	void Paint_Density_Field(const Grid & density, int screenid);
	void Paint_Velocity_Field(const Grid & U, const Grid & V, int screenid);
	void Paint(void);
	void Step();
	//void Apply_Particle_Force(const Particle & p);
	Float Get_U(Float x, Float y);
	Float Get_V(Float x, Float y);
	void Apply_Boundary_Condition(void);
	int U_Linear_ID(int i, int j);
	int V_Linear_ID(int i, int j);
	int P_Linear_ID(int i, int j);
	void Load_Pressure_Mat(int & eqid, vector<Tri>& A, Eigen::VectorXd & b);
	void Fill_Back_Solution(const VectorXd & b);
	Grid Advect_U(void);
	Grid Advect_V(void);
	Grid Advect(const Grid & f);
	void Apply_External_Force(void);
	VectorXd Project_Pressure(void);
	void Fill_Back_Pressure(const VectorXd & p_vec);
	void Pressure_Update(void);
	void Velocity_Step(void);
	void Apply_Particles();
	void Step_Fluid(void);
};