#include "solver.h"

Dye::Dye(int _r, int _g, int _b) {
	rgb[0] = _r;
	rgb[1] = _g;
	rgb[2] = _b;
	dens.resize(MESHH, MESHW);
	dens.setZero();
	src.resize(dens);
}

void GridCoor_to_ClipCoor(Float &i, Float &j, int H, int W, int screenid) {
	Float y = H - i + screenid * H;
	y = y * 2 / (H*TOTAL_SCREEN) - 1.0f;
	Float x = j * 2 / W - 1.0f;
	i = x;
	j = y;
}

void Show2Mesh(Float &i, Float &j) {
	i /= (SHOWH + 0.0) / (MESHH + 0.0);
	j /= (SHOWW + 0.0) / (MESHW + 0.0);
}

void Mesh2Clip(Float &i, Float &j, int screenid) {
	GridCoor_to_ClipCoor(i, j, MESHH, MESHW, screenid);
}

void Show2Clip(Float &i, Float &j, int screenid) {
	GridCoor_to_ClipCoor(i, j, SHOWH, SHOWW, screenid);
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
	CU.resize(U);
	CV.resize(V);
	P.resize(MESHH, MESHW);
	P.setZero();
	div.resize(MESHH, MESHW);
	div.setZero();
}

void Solver::Draw_Colors(int screenid) {
	glBegin(GL_POINTS);
	for (int is = 0; is < SHOWH; is++) {
		for (int js = 0; js < SHOWW; js++) {
			Float xm = is, ym = js;
			Show2Mesh(xm, ym);
			Float x = is, y = js;
			Show2Clip(x, y, screenid);
			Float rgb[3] = { 0,0,0 };
			for (int c = 0; c < colors.size(); c++) {
				Truncate_Position(colors[c].dens, xm, ym);
				Float ds = Interpolate(colors[c].dens, xm, ym);
				for (int d = 0; d < 3; d++) {
					rgb[d] += colors[c].rgb[d] * ds;
				}
			}
			glColor3f(rgb[0], rgb[1], rgb[2]);
			//glVertex2i(js, is + (TOTAL_SCREEN - screenid)*SHOWH);
			glVertex2f(x, y);
		}
	}
	glEnd();
	glFlush();
}

void Solver::Draw_Velocity_Field(const Grid &U, const Grid &V, int screenid) {
	Float vel_scale = DT;
	int stride = 1;
	//cout << U << endl;
	Assert(U.rows() == MESHH && U.cols() == MESHW, "paint U size unmatch");
	Assert(V.rows() == MESHH && V.cols() == MESHW, "paint V size unmatch");
	for (int i = 0; i < MESHH; i += stride) {
		for (int j = 0; j < MESHW; j += stride) {
			Float x0 = i, y0 = j;
			Float cx = U(i, j), cy = V(i, j);
			Float x1 = x0 + cx * vel_scale, y1 = y0 + cy * vel_scale;
			Mesh2Clip(x0, y0, screenid);
			Mesh2Clip(x1, y1, screenid);
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

void Solver::Draw(void) {
	Draw_Velocity_Field(U, V, 0);
	Draw_Colors(1);
}

void Solver::Step() {
	Step_Fluid();
}


