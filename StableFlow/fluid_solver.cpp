#include "fluid_solver.h"

Float Phi(Float r) {
	if (r > 2 || r < -2) return 0;
	else {
		cout << (1 + cos(PI*r / 2)) / 4 << endl;
		return (1 + cos(PI*r / 2)) / 4;
	}
}

Float Phi_XY(Float x, Float y) {
	return Phi(x / H)*Phi(y / H) / (H*H);
}


Float Solver::Get_U(Float x, Float y) {
	return Interpolate(U, x + 0.5, y);
}

Float Solver::Get_V(Float x, Float y) {
	return Interpolate(V, x, y + 0.5);
}


Grid Solver::Advect_U(void) {
	int n = U.rows(), m = U.cols();
	Grid u1(n, m);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			Float x = i - 0.5, y = j;
			Float xmid = x - 0.5*DT*Get_U(i, j);
			Float ymid = y - 0.5*DT*Get_V(i, j);
			Float x0 = i - DT * Get_U(xmid, ymid);
			Float y0 = i - DT * Get_V(xmid, ymid);
			u1(i, j) = Get_U(x0, y0);
		}
	}
	return u1;
}

Grid Solver::Advect_V(void) {
	int n = V.rows(), m = V.cols();
	Grid v1(n, m);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			Float x = i, y = j - 0.5;
			Float xmid = x - 0.5*DT*Get_U(i, j);
			Float ymid = y - 0.5*DT*Get_V(i, j);
			Float x0 = i - DT * Get_U(xmid, ymid);
			Float y0 = i - DT * Get_V(xmid, ymid);
			v1(i, j) = Get_V(x0, y0);
		}
	}
	return v1;
}

int Solver::U_Linear_ID(int i, int j) {
	Assert(0 <= i && i < U.rows(), "U_Linear_ID i out of range");
	Assert(0 <= j && j < U.cols(), "U_Linear_ID j out of range");
	return i * U.cols() + j;
}

int Solver::V_Linear_ID(int i, int j) {
	Assert(0 <= i && i < V.rows(), "V_Linear_ID i out of range");
	Assert(0 <= j && j < V.cols(), "V_Linear_ID j out of range");
	return U.size() + i * V.cols() + j;
}

int Solver::P_Linear_ID(int i, int j) {
	Assert(0 <= i && i < P.rows(), "P_Linear_ID i out of range");
	Assert(0 <= j && j < P.cols(), "P_Linear_ID j out of range");
	return i * P.cols() + j;
	//return U.size() + V.size() + i * P.cols() + j;
}

void Solver::Load_Pressure_Mat(int &eqid, vector<Tri> &A, Eigen::VectorXd &b) {
	int n = U.rows(), m = U.cols();
	Assert(V.rows() == n, "U and V rows not match");
	Assert(V.cols() == m, "U and V cols not match");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			Float rhs = 0;
			//p(i,j)
			A.push_back(Tri(eqid, P_Linear_ID(i, j), 2 * DT / (rho*DX*DX) + 2 * DT / (rho*DY*DY)));
			//p(i+1,j)
			if (i + 1 < n) A.push_back(Tri(eqid, P_Linear_ID(i + 1, j), -DT / (rho*DX*DX)));
			else A.push_back(Tri(eqid, P_Linear_ID(i + 1 - n, j), -DT / (rho*DX*DX)));
			//p(i-1,j)
			if (i - 1 >= 0) A.push_back(Tri(eqid, P_Linear_ID(i - 1, j), -DT / (rho*DX*DX)));
			else A.push_back(Tri(eqid, P_Linear_ID(i - 1 + n, j), -DT / (rho*DX*DX)));
			//p(i,j+1)
			if (j + 1 < m) A.push_back(Tri(eqid, P_Linear_ID(i, j + 1), -DT / (rho*DY*DY)));
			else A.push_back(Tri(eqid, P_Linear_ID(i, j + 1 - m), -DT / (rho*DY*DY)));
			//p(i,j-1)
			if (j - 1 >= 0) A.push_back(Tri(eqid, P_Linear_ID(i, j - 1), -DT / (rho*DY*DY)));
			else A.push_back(Tri(eqid, P_Linear_ID(i, j - 1 + m), -DT / (rho*DY*DY)));
			//u(i+1,j)
			if (i + 1 < n) rhs -= U(i + 1, j) / DX;
			else rhs -= U(i + 1 - n, j) / DX;
			//u(i, j)
			rhs += U(i, j) / DX;
			//v(i,j+1)
			if (j + 1 < m) rhs -= V(i, j + 1) / DY;
			else rhs -= V(i, j + 1 - m) / DY;
			//v(i,j)
			rhs += V(i, j) / DY;
			b(eqid) = rhs;
			eqid++;
		}
	}
}



void Solver::Fill_Back_Solution(const VectorXd &x) {
	int n = U.rows(), m = U.cols();
	Assert(V.rows() == n && V.cols() == m, "V,shape doesn't match");
	Assert(P.rows() == n && P.cols() == m, "P.shape doesn't match");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			U(i, j) = x(U_Linear_ID(i, j));
			V(i, j) = x(V_Linear_ID(i, j));
			P(i, j) = x(P_Linear_ID(i, j));
		}
	}
}


Grid Solver::Advect(const Grid &f) { // advect based on U and V
	int n = f.rows(), m = f.cols();
	Grid f1(n, m);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			Float xmid = i - 0.5*DT*U(i, j);
			Float ymid = j - 0.5*DT*V(i, j);
			Float x0 = i - DT * Interpolate(U, xmid, ymid);
			Float y0 = i - DT * Interpolate(V, xmid, ymid);
			f1(i, j) = Interpolate(f, x0, y0);
		}
	}
	return f1;
}

void Solver::Apply_External_Force(void) {
	for (int i = 0; i < U.rows(); i++) {
		for (int j = 0; j < U.cols(); j++) {
			U(i, j) += FU(i, j)*DT;
		}
	}
	for (int i = 0; i < V.rows(); i++) {
		for (int j = 0; j < V.cols(); j++) {
			V(i, j) += FV(i, j)*DT;
		}
	}
}

VectorXd Solver::Project_Pressure(void) {
	int eqid = 0;
	int num_eqs = U.size();
	vector<Tri> A_vec;
	VectorXd b(num_eqs);
	Load_Pressure_Mat(eqid, A_vec, b);
	SpMat A(num_eqs, num_eqs);
	A.setFromTriplets(A_vec.begin(), A_vec.end());
	//SpMat c(num_eqs, num_eqs);
	//c= A.transpose();
	//cout << "project matrix: " << endl << c-A<< endl;
	
	/*SparseLU<SparseMatrix<Float>, COLAMDOrdering<int> >   solver;
	solver.analyzePattern(A);
	solver.factorize(A);
	VectorXd X = solver.solve(b);*/
	
	ConjugateGradient<SparseMatrix<Float>, Lower | Upper> pcg;
	pcg.compute(A);
	VectorXd X = pcg.solve(b);
	//cout << "A=:" << endl << A << endl;
	//cout << "b=" << endl << b << endl;
	//cout << A << " " << b << endl;
	//VectorXd r = A * X - b;
	//cout << r.cwiseAbs().maxCoeff() << endl;
	return X;
}

void Solver::Fill_Back_Pressure(const VectorXd &p_vec) {
	int n = P.rows(), m = P.cols();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			P(i, j) = p_vec(P_Linear_ID(i, j));
		}
	}
}

void Solver::Pressure_Update(void) {
	int n = U.rows(), m = U.cols();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			//p(i,j)->u
			U(i, j) -= DT / (rho*DX)*P(i, j);
			//p(i-1,j)->u
			if (i - 1 >= 0) U(i, j) += DT / (rho*DX)*P(i - 1, j);
			else U(i, j) += DT / (rho*DX)*P(i - 1 + n, j);
			//p(i,j)->v
			V(i, j) -= DT / (rho*DY)*P(i, j);
			//p(i,j-1)->v
			if (j - 1 >= 0) V(i, j) += DT / (rho*DY)*P(i, j - 1);
			else V(i, j) += DT / (rho*DY)*P(i, j - 1 + m);
		}
	}
}

void Solver::Velocity_Step(void) {
	Apply_External_Force();
	U1 = Advect_U();
	V1 = Advect_V();

}

void Solver::Step_Fluid(void) {
	Grid U0 = U;
	Grid V0 = V;
	step_time++;
	printf("step fluid\n");
	U1 = Advect_U();
	V1 = Advect_V();
	//cout << "after advect V: " << V << endl;
	U = U1;
	V = V1;
	Apply_External_Force();
	//cout << "after apply external force V:" << V << endl;
	VectorXd p_vec = Project_Pressure();
	Fill_Back_Pressure(p_vec);
	//cout << "projected pressure: " << P << endl;
	Pressure_Update();
	//cout << (Grid_D(U, X, FORWARD) + Grid_D(V, Y, FORWARD)).abs().maxCoeff();
	//cout << "after projection V:" << V << endl;

	//cout << U.topRows(1) << endl;
	//cout << (Grid_D(U, X, CENTER) + Grid_D(V, Y, CENTER)).topRows(1) << endl;
	Grid u_res = rho * ((U - U0) / DT + Grid_S(U0, V0, U0)) - Grid_D(P, X, BACKWARD)
		- miu * (Grid_D(Grid_D(U, X, BACKWARD), X, FORWARD) + Grid_D(Grid_D(U, Y, BACKWARD), Y, FORWARD))
		- FU;
	cout << u_res << endl;
}