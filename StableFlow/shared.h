#pragma once

#define OPENGL
#define OPENMP
#define NDEBUG
//#define DIFFUSION_ON

#include <iostream>
#include <cstdio>
#include <tuple>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <list>
#include <float.h>
#include <ctime>
#include <omp.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#ifdef OPENGL
#include <GL/glut.h>
#endif

//NOTE!!! GL goes first, THEN Freeimage.
#include "FreeImage.h"
#pragma comment(lib,"FreeImage.lib")

using namespace std;
using namespace Eigen;
typedef Eigen::ArrayXXd Grid;
typedef char BYT;
typedef double Float;

const Float PI = acos(-1.0);

// Parameters about physics
const int MESHW = 200, MESHH = 100;
const int FLUID_PTNUM = (MESHW - 2)*(MESHH - 2);
const Float WHRATIO = MESHW / MESHH;
const Float DX = 0.01;
const Float DY = 0.01;
const Float H = 0.01;
const Float DT = 0.1;
const Float rho = 1;
const Float visc = 0;
const Float diff = 0;

// Parameters about solver
const int JACOBI_STEP = 20;
const int MAX_PCG_STEP = 10000;
const Float MFPCG_TOLERANCE = 1e-15;


// enumerations
enum DIFFTYPE { FORWARD, BACKWARD, CENTER };
enum MASKTYPE { VALID, INVALID };
enum AXIS { X, Y, N }; //N=neutral

// Parameters of rendering
const int TOTAL_SCREEN = 2;
const Float SCALE = 2;
const int SHOWW = MESHW*SCALE, SHOWH = MESHH*SCALE;
const int SHOW_FPS = 20;
const int TOTAL_FRAME = SHOW_FPS * 30;

void Assert(bool condition, const char *s);


