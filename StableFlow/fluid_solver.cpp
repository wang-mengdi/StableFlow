#include "fluid_solver.h"


void Apply_Boundary_Condition(Grid &F, AXIS ax) {
	int n = F.rows(), m = F.cols();
	for (int i = 1; i < m - 1; i++) {
		int sgn = (ax == X) ? -1 : 1;
		F(0, i) = sgn * F(1, i);
		F(n - 1, i) = sgn * F(n - 2, i);
	}
	for (int i = 1; i < n - 1; i++) {
		int sgn = (ax == Y) ? -1 : 1;
		F(i, 0) = sgn * F(i, 1);
		F(i, m - 1) = sgn * F(i, m - 2);
	}
	F(0, 0) = -F(1, 1);
	F(0, m - 1) = -F(1, m - 2);
	F(n - 1, 0) = -F(n - 2, 1);
	F(n - 1, m - 1) = -F(n - 2, m - 2);
}

void Jacobi_Solve(Grid &x, const Grid &b, Float a, Float c, AXIS ax) {
	// solve poisson function: nabla(x)=b within the interior of grid x
	int n = x.rows(), m = x.cols();
	for (int tot = 0; tot < JACOBI_STEP; tot++) {
		Float mxdiff = 0;
		for (int i = 1; i < n - 1; i++) {
			for (int j = 1; j < m - 1; j++) {
				Float x1 = (b(i, j) + a * (x(i - 1, j) + x(i + 1, j) + x(i, j - 1) + x(i, j + 1))) / c;
				mxdiff = max(mxdiff, fabs(x(i, j) - x1));
				x(i, j) = x1;
			}
		}
		Apply_Boundary_Condition(x, ax);
	}
}

void Solver::Diffuse(AXIS ax, Grid &x, const Grid &b, Float diff) {
	Float a = DT * diff;
	Jacobi_Solve(x, b, a, 1 + 4 * a, ax);
}

Grid SL_Advect(AXIS ax, const Grid &q, const Grid &U, const Grid &V, DIFFTYPE dir) {
	//Semi Lagrangian Advect
	int n = q.rows(), m = q.cols();
	Grid q1;
	q1.resize(n, m);
	q1.setZero();
	for (int i = 1; i < n - 1; i++) {
		for (int j = 1; j < m - 1; j++) {
			Float x, y;
			if (dir == BACKWARD) {
				x = i - DT * U(i, j), y = j - DT * V(i, j);
			}
			else if (dir == FORWARD) {
				x = i + DT * U(i, j), y = j + DT * V(i, j);
			}
			x = max(x, 0.5);
			x = min(x, n - 1.5);
			y = max(y, 0.5);
			y = min(y, m - 1.5);
			q1(i, j) = Interpolate(q, x, y);
		}
	}
	Apply_Boundary_Condition(q1, ax);
	return q1;
}

Grid MacCormack(AXIS ax, const Grid &q, const Grid &U, const Grid &V) {
	Grid q_est = SL_Advect(ax, q, U, V, BACKWARD);
	Grid q_bck = SL_Advect(ax, q_est, U, V, FORWARD);
	Grid err = q - q_bck;
	return q_est + 0.5*err;
}

int Solve_Field_ID(const Grid &A, int i, int j) {
	int n = A.rows(), m = A.cols();
	Assert(1 <= i && i < n - 1, "i out of bound");
	Assert(1 <= j && j < m - 1, "j out of bound");
	i--, j--;
	return i * (m - 2) + j;
}

void Load_Pressure_System(SpMat &A, Eigen::VectorXd &b, const Grid &div) {
	const int dx[4] = { 1,-1,0,0 }, dy[4] = { 0,0,1,-1 };
	int n = div.rows(), m = div.cols();
	int numeq = (n - 2)*(m - 2);
	Assert(A.rows() == numeq && A.cols() == numeq, "sparse matrix size not match");
	vector<Tri> v;
	for (int i = 1; i < n - 1; i++) {
		for (int j = 1; j < m - 1; j++) {
			Float ctr = -4;
			int eid = Solve_Field_ID(div, i, j);
			for (int d = 0; d < 4; d++) {
				int i1 = i + dx[d], j1 = j + dy[d];
				if (i1 == 0 || i1 == n - 1 || j1 == 0 || j1 == m - 1) {//boundary, it's the same as (i,j)
					ctr += 1;
				}
				else {
					int id1 = Solve_Field_ID(div, i1, j1);
					v.push_back(Tri(eid, id1, 1));
				}
			}
			v.push_back(Tri(eid, eid, ctr));
			b(eid) = div(i, j);
		}
	}
	A.setFromTriplets(v.begin(), v.end());
}

void Fill_Pressure_Solution(Grid &P, const VectorXd &x) {
	int n = P.rows(), m = P.cols();
	for (int i = 1; i < n - 1; i++) {
		for (int j = 1; j < m - 1; j++) {
			int eid = Solve_Field_ID(P, i, j);
			P(i, j) = x(eid);
		}
	}
	Apply_Boundary_Condition(P, N);
}

void Get_Div(const Grid &U, const Grid &V, Grid &div) {
	int n = div.rows(), m = div.cols();
	for (int i = 1; i < n - 1; i++) {
		for (int j = 1; j < m - 1; j++) {
			div(i, j) = 0.5f*(U(i + 1, j) - U(i - 1, j) + V(i, j + 1) - V(i, j - 1));
		}
	}
	Apply_Boundary_Condition(div, N);
}

Grid Laplace(const Grid &P) {
	int n = P.rows(), m = P.cols();
	Grid lp(n, m);
	for (int i = 1; i < n - 1; i++) {
		for (int j = 1; j < m - 1; j++) {
			lp(i, j) = P(i + 1, j) + P(i - 1, j) + P(i, j - 1) + P(i, j + 1) - 4 * P(i, j);
		}
	}
	Apply_Boundary_Condition(lp, N);
	return lp;
}

void Get_Pressure(Grid &P, const Grid &U, const Grid &V, Grid &div) {
	int n = P.rows(), m = P.cols();
	Get_Div(U, V, div);
	Apply_Boundary_Condition(div, N);

	int numeq = (n - 2)*(m - 2);
	SpMat A(numeq, numeq);
	VectorXd b(numeq);
	Load_Pressure_System(A, b, div);
	ConjugateGradient<SparseMatrix<Float>, Lower | Upper > pcg;
	pcg.compute(A);
	VectorXd X = pcg.solve(b);
	VectorXd r = A * X - b;
	//cout << "PCG residual norm: " << r.norm() << endl;
	Fill_Pressure_Solution(P, X);

	/*P.setZero();
	Jacobi_Solve(P, -div, 1, 4, N);*/

	Grid L = Laplace(P);
	//cout << "lap and grad diff norm: " << Grid_Norm(L - div) << endl;
	cout << "solved P 2norm: " << Grid_Norm(P) << endl;
}

void Update_Pressure(Grid &U, Grid &V, const Grid &P) {
	int n = P.rows(), m = P.cols();
	for (int i = 1; i < n - 1; i++) {
		for (int j = 1; j < m - 1; j++) {
			U(i, j) -= (P(i + 1, j) - P(i - 1, j)) / 2;
			V(i, j) -= (P(i, j + 1) - P(i, j - 1)) / 2;
		}
	}
	Apply_Boundary_Condition(U, X);
	Apply_Boundary_Condition(V, Y);
}

void Project(Grid &U, Grid &V, Grid &P, Grid &div) {
	Get_Pressure(P, U, V, div);
	Update_Pressure(U, V, P);
}

void Add_Source(Grid &f, const Grid &s) {
	int n = f.rows(), m = f.cols();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			f(i, j) = s(i, j);
		}
	}
}

void Solver::Velocity_Step(void) {
	Grid U0, V0;
	Apply_ConstMask(U, CU);
	Apply_ConstMask(V, CV);
#ifdef DIFFUSION_ON
	U0 = U;
	V0 = V;
	Diffuse(X, U, U0, visc);
	Diffuse(Y, V, V0, visc);
	Project(U, V, P, div);
#endif
	U0 = U, V0 = V;
	U = SL_Advect(X, U, U0, V0, BACKWARD);
	V = SL_Advect(Y, V, U0, V0, BACKWARD);
	//U = MacCormack(X, U, U0, V0);
	//V = MacCormack(Y, V, U0, V0);
	Project(U, V, P, div);
}

void Solver::Density_Step(Dye &D) {
	Apply_ConstMask(D.dens, D.src);
	//Grid D0 = D.dens;
#ifdef DIFFUSION_ON
	Diffuse(N, D.dens, D0, diff);
#endif
	//D0 = D.dens;
	//D.dens = MacCormack(N, D.dens, U, V);
	D.dens = SL_Advect(N, D.dens, U, V, BACKWARD);
}

void Solver::Step_Fluid(void) {
	step_time++;
	cout << "fluid step time " << step_time << endl;
	Velocity_Step();
	for (int i = 0; i < colors.size(); i++) {
		Density_Step(colors[i]);
	}
	//cout << "density: " << colors[0].dens << endl;
	Get_Div(U, V, div);
	cout << "step done, divergence norm:" << Grid_Norm(div) << endl;
}
