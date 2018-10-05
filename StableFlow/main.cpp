#include "shared.h"
#include "solver.h"
Solver S;

void Init(void) {
	Float r = 0.04;
	Float v = 30;
	Float h = 0.17;
	Float x1 = 0.5 - h, y1 = 0.1;
	Float x2 = 0.5 + h, y2 = 0.1;
	S.Init();
	Dye dred = Dye(1, 0, 0);
	dred.src.Set_Real_Circle(x1,y1, r, 1);
	S.colors.push_back(dred);
	Dye dblue = Dye(0, 0, 1);
	dblue.src.Set_Real_Circle(x2,y2, r, 1);
	S.colors.push_back(dblue);
	S.CV.Set_Real_Circle(x1, y1, r, v + 1);
	S.CU.Set_Real_Circle(x1, y1, r, v);
	S.CV.Set_Real_Circle(x2, y2, r, v);
	S.CU.Set_Real_Circle(x2, y2, r, -v);
	//S.CV.Set_Box(0.4, 0.6, 0.2, 0.3, 10);
	//S.CV.Set_Box(0.4, 0.6, 0.7, 0.8, -10);
}


void Process_Click(int button, int state, int y, int x) {
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
	}
}

void Plot(void) {
	glClear(GL_COLOR_BUFFER_BIT);
	S.Draw();
	glFlush();
	glutSwapBuffers();
}



void Step_Time(int v) {
	S.Step();
	Plot();
	glutTimerFunc(1000 / SHOW_FPS, Step_Time, 1);
}

int main(int argc, char *argv[]) {
	Init();
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);//double buffer
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(SHOWW, SHOWH*TOTAL_SCREEN);
	int glut_window = glutCreateWindow("Rasterization");
	glutMouseFunc(&Process_Click);
	glutDisplayFunc(&Plot);
	glutTimerFunc(1000 / SHOW_FPS, Step_Time, 1);
	glutMainLoop();
	return 0;
}