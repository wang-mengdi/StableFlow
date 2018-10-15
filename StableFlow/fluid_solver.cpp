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

void Advect(AXIS ax, Grid &q, const Grid &q0, const Grid &U, const Grid &V) {
	int n = q.rows(), m = q.cols();
	for (int i = 1; i < n - 1; i++) {
		for (int j = 1; j < m - 1; j++) {
			Float x = i - DT * U(i, j), y = j - DT * V(i, j);
			x = max(x, 0.5);
			x = min(x, n - 1.5);
			y = max(y, 0.5);
			y = min(y, m - 1.5);
			q(i, j) = Interpolate(q0, x, y);
		}
	}
	Apply_Boundary_Condition(q, ax);
}

void Get_Div(const Grid &U, const Grid &V, Grid &div) {
	int n = div.rows(), m = div.cols();
	for (int i = 1; i < n - 1; i++) {
		for (int j = 1; j < m - 1; j++) {
			div(i, j) = -0.5f*(U(i + 1, j) - U(i - 1, j) + V(i, j + 1) - V(i, j - 1));
		}
	}
}

void Get_Pressure(Grid &P, const Grid &U, const Grid &V, Grid &div) {
	int n = P.rows(), m = P.cols();
	Get_Div(U, V, div);
	Apply_Boundary_Condition(div, N);
	P.setZero();
	Jacobi_Solve(P, div, 1, 4, N);
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
	Advect(X, U, U0, U0, V0);
	Advect(Y, V, V0, U0, V0);
	Project(U, V, P, div);
}

void Solver::Density_Step(Dye &D) {
	Apply_ConstMask(D.dens, D.src);
	Grid D0 = D.dens;
	Diffuse(N, D.dens, D0, diff);
	D0 = D.dens;
	Advect(N, D.dens, D0, U, V);
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
	cout << "step done, max divergence:" << div.abs().maxCoeff() << endl;
}
