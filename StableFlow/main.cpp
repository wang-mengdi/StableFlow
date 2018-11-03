#include "shared.h"
#include "solver.h"
Solver S;

void Init(void) {
	FreeImage_Initialise(TRUE);
	Float r = 0.04;
	Float v = 30;
	Float h = 0.17;
	Float x1 = 0.5 - h, y1 = 0.1;
	Float x2 = 0.5 + h, y2 = 0.1;
	S.Init();
	Dye dred = Dye(1, 0, 0);
	dred.src.Set_Real_Circle(x1,y1, r, 1);
	S.colors.push_back(dred);
	S.CV.Set_Real_Circle(x1, y1, r, v + 1);
	S.CU.Set_Real_Circle(x1, y1, r, v);
	Dye dblue = Dye(0, 0, 1);
	dblue.src.Set_Real_Circle(x2,y2, r, 1);
	S.colors.push_back(dblue);
	S.CV.Set_Real_Circle(x2, y2, r, v);
	S.CU.Set_Real_Circle(x2, y2, r, -v);
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

void Save_Image(int step) {
	int width = SHOWW;
	int height = SHOWH * TOTAL_SCREEN;

	static BYTE *pixels = (BYTE*)malloc(width * height * 3);

	glReadPixels(0, 0, width, height, GL_BGR_EXT, GL_UNSIGNED_BYTE, pixels);

	FIBITMAP *image = FreeImage_ConvertFromRawBits(pixels, width, height, 3 * width, 24, FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK, FALSE);

	static char filename[256];
	sprintf_s(filename, 256, "./images/%08d.png", step);

	if (FreeImage_Save(FIF_PNG, image, filename, 0))
		printf("Successfully saved!\n");
	else
		printf("Failed saving!\n");

	FreeImage_Unload(image);

	//free(pixels);
}

void Step_Time(int v) {
	S.Step();
	Plot();
	Save_Image(S.step_time);
	glutTimerFunc(1000 / SHOW_FPS, Step_Time, 1);
	if (S.step_time >= TOTAL_FRAME) {
		exit(0);
	}
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