#include "solver.h"

void GridCoor_to_ClipCoor(Float &x, Float &y, int screenid) {
	y = MESHH - y + screenid * MESHH;
	x = x * 2 / MESHW - 1.0f;
	y = y * 2 / (MESHH*TOTAL_SCREEN) - 1.0f;
}

void ScreenCoor_to_ClipCoor(Float &x, Float &y, int screenid) {
	x /= SCALE;
	y /= SCALE;
	y = MESHH - y + screenid * MESHH;
	x = x * 2 / MESHW - 1.0f;
	y = y * 2 / (MESHH*TOTAL_SCREEN) - 1.0f;
}

void Solver::Init(void) {
	U.resize(MESHH, MESHW);
	U.setZero();
	V.resize(MESHH, MESHW);
	V.setZero();
	U1.resize(MESHH, MESHW);
	U1.setZero();
	V1.resize(MESHH, MESHW);
	V1.setZero();
	P.resize(MESHH, MESHW);
	P.setZero();
	FU.resize(MESHH, MESHW);
	FU.setZero();
	FV.resize(MESHH, MESHW);
	FV.setZero();
}

void Solver::Paint_Density_Field(const Grid &density, int screenid) {
	glBegin(GL_POINTS);
	Float mx = -1;
	for (int j = 0; j < MESHH; j++) {
		for (int k = 0; k < MESHW; k++) {
			mx = max(mx, fabs(density(j, k)));
		}
	}
	mx = max(mx, (Float)0.0005);
	cout << "max density: "<<mx << endl;
	for (int is = 0; is < MESHH*SCALE; is++) {
		for (int js = 0; js < MESHW*SCALE; js++) {
			int i = (is + 0.0) / SCALE;
			int j = (js + 0.0) / SCALE;
			Float x = is, y = js;
			ScreenCoor_to_ClipCoor(x, y, screenid);
			Float d = fabs(density(i, j) / mx);
			//if (d != 0) printf("%f ", d);
			glColor3f(d, d, d);
			glVertex2f(x, y);
		}
	}
	glEnd();
	glFlush();
}

void Solver::Paint_Velocity_Field(const Grid &U, const Grid &V, int screenid) {
	Float vel_scale = 10;
	int stride = 1;
	//cout << U << endl;
	Assert(U.rows() == MESHH && U.cols() == MESHW, "paint U size unmatch");
	Assert(V.rows() == MESHH && V.cols() == MESHW, "paint V size unmatch");
	for (int i = 0; i < MESHH; i += stride) {
		for (int j = 0; j < MESHW; j += stride) {
			Float x0 = j, y0 = i;
			Float cx = V(i, j), cy = U(i, j);
			Float x1 = x0 + cx * vel_scale, y1 = y0 - cy * vel_scale;
			GridCoor_to_ClipCoor(x0, y0, screenid);
			GridCoor_to_ClipCoor(x1, y1, screenid);
			glColor3f(0.0, 1.0, 0.0);
			glBegin(GL_POINTS);
			glVertex2f(x0, y0);
			glEnd();
			glColor3f(1.0, 1.0, 1.0);
			glBegin(GL_LINES);
			glVertex2f(x0, y0);
			glVertex2f(x1, y1);
			glEnd();
		}
	}
}

void Solver::Paint(void) {
	Paint_Velocity_Field(U, V, 0);
	Paint_Density_Field(P, 1);
	Paint_Velocity_Field(U1, V1, 2);
}

void Solver::Step() {
	//Apply_Particles(particles);
	Step_Fluid();
}