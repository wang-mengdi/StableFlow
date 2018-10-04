#include "shared.h"
#include "solver.h"
//#include "fiber.h"

Solver S;

//vector<Particle> Particles;

void Init(void) {
	S.Init();
	S.FV.block(4, 2, 1, 1) = 1;
	S.FV.block(4, 6, 1, 1) = -1;
}


void Process_Click(int button, int state, int y, int x) {
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
		printf("add blast\n");
		//balls.push_back(BlastBall(x, y));
	}
}

void Plot(void) {
	glClear(GL_COLOR_BUFFER_BIT);
	//S.Paint();
	glFlush();
	glutSwapBuffers();
}



void Step_Time(int v) {
	//S.Step(Particles);
	//Plot();
	glutTimerFunc(1000 / SHOW_FPS, Step_Time, 1);
}

int main(int argc, char *argv[]) {
	//Init();
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