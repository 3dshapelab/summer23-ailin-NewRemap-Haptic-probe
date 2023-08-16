// This script builds the visual stimulus for haptic remapping experiment for summer 2023
// to change the display distance, change it from the variable declaration
// to save time, we can set bool resetScreen_betweenRuns to false

#include <cstdlib>
#include <cmath>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <random>

/**** BOOOST MULTITHREADED LIBRARY *********/
#include <boost/thread/thread.hpp>
#include <boost/asio.hpp>	//include asio in order to avoid the "winsock already declared problem"


#ifdef _WIN32
#include <windows.h>
#include <gl\gl.h>            // Header File For The OpenGL32 Library
#include <gl\glu.h>            // Header File For The GLu32 Library
#include "glut.h"            // Header File For The GLu32 Library
#include <MMSystem.h>
#endif

/************ INCLUDE CNCSVISION LIBRARY HEADERS ****************/
//#include "Optotrak.h"
//#include "Optotrak2.h"
#include "Marker.h"
#include "Mathcommon.h"
#include "GLUtils.h"
#include "VRCamera.h"
#include "CoordinatesExtractor.h"
#include "StimulusDrawer.h"
#include "GLText.h"
#include "Util.h"
#include "SOIL.h"
#define BROWN 
#ifdef BROWN
#include "BrownMotorFunctions.h"
#else
#include "RoveretoMotorFunctions.h"
#endif
/********* NAMESPACE DIRECTIVES ************************/
using namespace std;
using namespace mathcommon;
using namespace Eigen;
using namespace util;
using namespace BrownMotorFunctions;

/********* #DEFINE DIRECTIVES **************************/
#include "Calibration_017B.h"
const float DEG2RAD = M_PI / 180;
Screen screen;
/********* VISUALIZATION VARIABLES *****************/
static const bool gameMode = true;
bool stimulus_built = false;
bool visibleInfo = true;

/********* VARIABLES OBJECTS  **********************/
VRCamera cam;
//Optotrak2 *optotrak;
CoordinatesExtractor headEyeCoords;


/********** EYES AND MARKERS **********************/
Vector3d eyeLeft, eyeRight, eyeCyclopean;
vector <Marker> markers;
double interoculardistance = 60;

/********** STIMULUS SHAPE ***************/
// display distance and size
double display_distance = -400;
double dist_toEye = -display_distance;
double visual_angle = 7.4; // stim diangonal size

//height and width
double stimulus_height = 70; //tan((DEG2RAD * visual_angle)/2) * 2 * (abs(display_distance));
double stimulus_width = 70; //ratio_bgwidth_height * stimulus_height;
double stimulus_visiblewidth = 70; //ratio_visiblewidth_height * stimulus_height;
double ratio_width_height = 1.2;//1.3;
double ratio_width_height_long = 1.4;
double ratio_visiblewidth_height = 1.15;//1.1;

// depths of visual stimuli
double depth_text = 32;
double depth_disp = 32;
double depth_inc = 2;

/********** STIMULUS VERTICES ***************/
struct Vec2 {
	float x, y;
};

struct CurveYLMap {
	std::vector<double> y_vec;
	std::vector<double> l_vec;
	double curve_depth;
	double curve_height;
	double step_size;
};

struct CurvePtsData {
	std::vector<double> y_vec;
	std::vector<double> z_vec;
	std::vector<double> l_vec;
	double curve_depth;
	double curve_height;
	double step_size;
};


struct TextureDotsData {
	std::vector<Vec2> dot_center_vec;
	Vec2 TexMapSize;
	float Radius;
	float margin_y;
};

struct VerticesData {
	std::vector<GLfloat> vertices_vec;
	std::vector<GLfloat> colors_vec;
	std::vector<GLfloat> light_normals_vec;
	std::vector<GLfloat> textuv_vec;
	std::vector<GLuint> indices_draw_triangle_vec;

};

struct ContourData {
	std::vector<Vector3f> vert_Rcontour;
	std::vector<Vector3f> vert_Lcontour;
};

VerticesData my_verts;
ContourData my_contour_data;

/********** RAMDOM DOTS SURF ***************/
// dot distribution
int dot_num_per_col = 18;
int dot_number;
double dot_jitter_max_scale = 1.30;

// for ot size
double visAngle_dot = 0.13;
double ratio_dotSize_distance = tan(DEG2RAD * visAngle_dot / 2.0);
int dot_sides = 16; // Bigger dot: 24; smaller dot: 8

// for random dot surface
float bgSurface_color = 0.8;
int nr_points_height_bgSurface = 76;

std::vector<Vector3d> dot_container;

/********** OPTIONS ***************/
bool resetScreen_betweenRuns = false;


/********** FUNCTION PROTOTYPES *****/
void beepOk(int tone);
void drawGLScene();
void handleKeypress(unsigned char key, int x, int y);
void handleResize(int w, int h);
void initProjectionScreen(double _focalDist, const Affine3d& _transformation = Affine3d::Identity(), bool synchronous = true);
void update(int value);
void idle();
void initMotors();
void initGLVariables();
void initVariables();
void initStreams();
void initRendering();
void initStimulus();
void drawInfo();
void drawStimulus();


double getZ(double shapeHeight, double shapeDepth, double Y);
void generateRandomDots(double shapeWidth, double shapeHeight, double shapeDepth, int dotNumPerRow, int dotNumPerCol, double dotJitterMax_Scale, vector<Vector3d>& dotContainer);
void buildContour_wavy(double ContourWidth, double shapeHeight, double shapeDepth, ContourData& new_contours_vert);
void buildRandomDotSurface(double shapeWidth, double shapeHeight, double shapeDepth, double contourPanelSeparation, VerticesData& vertices_data, vector<Vector3d>& dotContainer, ContourData& contours_vert);
void drawRandomDots(const vector<Vector3d>& dotContainer);
void drawRandomDotSurface(double distShapeToEye, const VerticesData& vertices_data, const vector<Vector3d>& dotContainer, const ContourData& contours_vert);
void drawContours(const ContourData& contours_vert);



double getZ(double shapeHeight, double shapeDepth, double Y) {

	double Z;

	Z = shapeDepth * cos(M_PI * Y / shapeHeight);

	return (Z);
}


void initStimulus() {

	stimulus_built = false;
	dist_toEye = -(display_distance - depth_disp);
	stimulus_height = tan((DEG2RAD * visual_angle) / 2) * 2 * dist_toEye;
	stimulus_visiblewidth = ratio_visiblewidth_height * stimulus_height;
	stimulus_width = ratio_width_height_long * stimulus_height;
	buildRandomDotSurface(stimulus_width, stimulus_height, depth_disp, stimulus_visiblewidth, my_verts, dot_container, my_contour_data);
	stimulus_built = true;
}

void buildContour_wavy(double ContourWidth, double shapeHeight, double shapeDepth, ContourData& new_contours_vert) {

	new_contours_vert = {};

	int nr_seg = 1, nr_wave_perSeg = 6, nr_unit_perWave = 8, nr_step_unit = 10; // number of segments need to be an even number
	int del_unit_prev = 0, del_unit_next = 0, nr_unit_current = nr_unit_perWave;

	double step_y = shapeHeight / (nr_seg * nr_wave_perSeg * nr_unit_perWave * nr_step_unit);

	int rand_portions = 20; //how many portions we want - will creat (1 + rand_portions) different values
	double rand_min = .5, rand_max = 2;
	double seg_x_off_L, seg_x_off_R;

	double x_L = -ContourWidth / 2, x_R = ContourWidth / 2; // 
	double y = -shapeHeight / 2, y0 = -shapeHeight / 2;
	double tempy_height = ((double)(1.0 / 6.0)) * (shapeHeight / nr_seg);
	double dummy_sign = -1;

	double z = 0;

	new_contours_vert.vert_Lcontour.push_back(Vector3d(x_L, y, z));
	new_contours_vert.vert_Rcontour.push_back(Vector3d(x_R, y, z));


	// Left
	for (int i_wave = 0; i_wave < nr_wave_perSeg; i_wave++) {

		if (i_wave < (nr_wave_perSeg - 1)) {
			del_unit_next = rand() % 5 - 2;
		}
		else {
			del_unit_next = 0;
		}

		nr_unit_current = nr_unit_perWave - del_unit_prev + del_unit_next;
		tempy_height = nr_unit_current * nr_step_unit * step_y;


		if (i_wave % 2 == 0) {

			dummy_sign = -1 * dummy_sign;
			seg_x_off_L = rand_min + (rand_max - rand_min) * (rand() % (rand_portions + 1)) / ((double)(rand_portions));
			//seg_x_off_R = rand_min + (rand_max - rand_min) * ( rand() % (rand_portions + 1) )/((double)(rand_portions));

			for (int i_step = 0; i_step < nr_unit_current * nr_step_unit; i_step++) {

				y = y + step_y;
				x_L = -ContourWidth / 2 - dummy_sign * seg_x_off_L * sin(0.5 * M_PI * (y - y0) / tempy_height);
				//x_R = ContourWidth / 2 + dummy_sign * seg_x_off_R * sin(0.5 * M_PI * (y - y0) / tempy_height);
				z = getZ(shapeHeight, shapeDepth, y);

				new_contours_vert.vert_Lcontour.push_back(Vector3d(x_L, y, z));

			}

		}
		else {

			for (int i_step = 0; i_step < nr_unit_current * nr_step_unit; i_step++) {

				y = y + step_y;
				x_L = -ContourWidth / 2 - dummy_sign * seg_x_off_L * cos(0.5 * M_PI * (y - y0) / tempy_height);
				//x_R = ContourWidth / 2 + dummy_sign * seg_x_off_R * cos(0.5 * M_PI * (y - y0) / tempy_height);	
				z = getZ(shapeHeight, shapeDepth, y);

				new_contours_vert.vert_Lcontour.push_back(Vector3d(x_L, y, z));

			}

		}

		y0 = y;
		del_unit_prev = del_unit_next;

	}

	// Right
	del_unit_prev = 0;
	y = -shapeHeight / 2, y0 = -shapeHeight / 2;
	dummy_sign = -1;
	for (int i_wave = 0; i_wave < nr_wave_perSeg; i_wave++) {

		if (i_wave < (nr_wave_perSeg - 1)) {
			del_unit_next = rand() % 5 - 2;
		}
		else {
			del_unit_next = 0;
		}

		nr_unit_current = nr_unit_perWave - del_unit_prev + del_unit_next;
		tempy_height = nr_unit_current * nr_step_unit * step_y;

		if (i_wave % 2 == 0) {

			dummy_sign = -1 * dummy_sign;
			//seg_x_off_L = rand_min + (rand_max - rand_min) * ( rand() % (rand_portions + 1) )/((double)(rand_portions));
			seg_x_off_R = rand_min + (rand_max - rand_min) * (rand() % (rand_portions + 1)) / ((double)(rand_portions));

			for (int i_step = 0; i_step < nr_unit_current * nr_step_unit; i_step++) {

				y = y + step_y;
				//x_L = -ContourWidth / 2 - dummy_sign * seg_x_off_L * sin(0.5 * M_PI * (y - y0) / tempy_height);
				x_R = ContourWidth / 2 + dummy_sign * seg_x_off_R * sin(0.5 * M_PI * (y - y0) / tempy_height);
				z = getZ(shapeHeight, shapeDepth, y);

				new_contours_vert.vert_Rcontour.push_back(Vector3d(x_R, y, z));

			}

		}
		else {

			for (int i_step = 0; i_step < nr_unit_current * nr_step_unit; i_step++) {

				y = y + step_y;
				//x_L = -ContourWidth / 2 - dummy_sign * seg_x_off_L * cos(0.5 * M_PI * (y - y0) / tempy_height);
				x_R = ContourWidth / 2 + dummy_sign * seg_x_off_R * cos(0.5 * M_PI * (y - y0) / tempy_height);
				z = getZ(shapeHeight, shapeDepth, y);

				new_contours_vert.vert_Rcontour.push_back(Vector3d(x_R, y, z));
			}

		}

		y0 = y;
		del_unit_prev = del_unit_next;

	}

}

void drawRandomDotSurface(double distShapeToEye, const VerticesData& vertices_data, const vector<Vector3d>& dotContainer, const ContourData& contours_vert) {

	glPushMatrix();

	glLoadIdentity();
	glTranslated(0, 0, -distShapeToEye);

	// activate and specify pointer to vertex array
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	glVertexPointer(3, GL_FLOAT, 0, &vertices_data.vertices_vec[0]);
	glColorPointer(3, GL_FLOAT, 0, &vertices_data.colors_vec[0]);
	glDrawElements(GL_TRIANGLES, vertices_data.indices_draw_triangle_vec.size(), GL_UNSIGNED_INT, &vertices_data.indices_draw_triangle_vec[0]);

	// deactivate vertex arrays after drawing
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);

	drawRandomDots(dotContainer);
	drawContours(contours_vert);

	glPopMatrix();
}


void drawContours(const ContourData& contours_vert) {

	int n;
	float panel_width = 40;
	float panel_height_extra = 20;

	glTranslated(0, 0, 2);
	n = int(contours_vert.vert_Lcontour.size());
	glColor3f(0.0f, 0.0f, 0.0f);
	if (n > 0) {

		// Right panels
		glBegin(GL_QUAD_STRIP);

		glVertex3f(contours_vert.vert_Rcontour.at(0)[0] + panel_width, contours_vert.vert_Rcontour.at(0)[1] - panel_height_extra, contours_vert.vert_Rcontour.at(0)[2]); //0
		glVertex3f(contours_vert.vert_Rcontour.at(0)[0], contours_vert.vert_Rcontour.at(0)[1] - panel_height_extra, contours_vert.vert_Rcontour.at(0)[2]); //1

		for (int i = 0; i < n; i++)
		{
			glVertex3f(contours_vert.vert_Rcontour.at(i)[0] + panel_width, contours_vert.vert_Rcontour.at(i)[1], contours_vert.vert_Rcontour.at(i)[2]); //0
			glVertex3f(contours_vert.vert_Rcontour.at(i)[0], contours_vert.vert_Rcontour.at(i)[1], contours_vert.vert_Rcontour.at(i)[2]); //1

		}

		glVertex3f(contours_vert.vert_Rcontour.at(n - 1)[0] + panel_width, contours_vert.vert_Rcontour.at(n - 1)[1] + panel_height_extra, contours_vert.vert_Rcontour.at(n - 1)[2]); //0
		glVertex3f(contours_vert.vert_Rcontour.at(n - 1)[0], contours_vert.vert_Rcontour.at(n - 1)[1] + panel_height_extra, contours_vert.vert_Rcontour.at(n - 1)[2]); //1

		glEnd();

		// Left panels
		glBegin(GL_QUAD_STRIP);

		glVertex3f(contours_vert.vert_Lcontour.at(0)[0], contours_vert.vert_Lcontour.at(0)[1] - panel_height_extra, contours_vert.vert_Lcontour.at(0)[2]); //0
		glVertex3f(contours_vert.vert_Lcontour.at(0)[0] - panel_width, contours_vert.vert_Lcontour.at(0)[1] - panel_height_extra, contours_vert.vert_Lcontour.at(0)[2]); //1

		for (int i = 0; i < n; i++)
		{
			glVertex3f(contours_vert.vert_Lcontour.at(i)[0], contours_vert.vert_Lcontour.at(i)[1], contours_vert.vert_Lcontour.at(i)[2]); //0
			glVertex3f(contours_vert.vert_Lcontour.at(i)[0] - panel_width, contours_vert.vert_Lcontour.at(i)[1], contours_vert.vert_Lcontour.at(i)[2]); //1

		}

		glVertex3f(contours_vert.vert_Lcontour.at(n - 1)[0], contours_vert.vert_Lcontour.at(n - 1)[1] + panel_height_extra, contours_vert.vert_Lcontour.at(n - 1)[2]); //0
		glVertex3f(contours_vert.vert_Lcontour.at(n - 1)[0] - panel_width, contours_vert.vert_Lcontour.at(n - 1)[1] + panel_height_extra, contours_vert.vert_Lcontour.at(n - 1)[2]); //1

		glEnd();
	}

}


void buildRandomDotSurface(double shapeWidth, double shapeHeight, double shapeDepth, double contourPanelSeparation, VerticesData& vertices_data, vector<Vector3d>& dotContainer, ContourData& contours_vert) {

	vertices_data = {};

	int nr_pts_height = nr_points_height_bgSurface;
	int nr_pts_width = nr_pts_height * (shapeWidth / shapeHeight);

	double step_size_width = shapeWidth / (double)(nr_pts_width - 1);
	double step_size_height = shapeHeight / (double)(nr_pts_height - 1);

	GLuint i_ind = 0;

	double x, y, z;
	//double normal_x, normal_y, normal_z;
	for (int jj = 0; jj < nr_pts_height; jj++) {

		y = -shapeHeight / 2 + jj * step_size_height;
		z = getZ(shapeHeight, shapeDepth, y);


		for (int ii = 0; ii < nr_pts_width; ii++) {
			x = -shapeWidth / 2 + ii * step_size_width;

			vertices_data.vertices_vec.push_back(x);
			vertices_data.vertices_vec.push_back(y);
			vertices_data.vertices_vec.push_back(z);

			vertices_data.colors_vec.push_back(bgSurface_color);
			vertices_data.colors_vec.push_back(0);
			vertices_data.colors_vec.push_back(0);

			if (ii < nr_pts_width - 1 && jj < nr_pts_height - 1) {

				// using vector
				vertices_data.indices_draw_triangle_vec.push_back(i_ind);
				vertices_data.indices_draw_triangle_vec.push_back(i_ind + 1);
				vertices_data.indices_draw_triangle_vec.push_back(i_ind + nr_pts_width);

				vertices_data.indices_draw_triangle_vec.push_back(i_ind + nr_pts_width);
				vertices_data.indices_draw_triangle_vec.push_back(i_ind + 1);
				vertices_data.indices_draw_triangle_vec.push_back(i_ind + nr_pts_width + 1);

			}

			i_ind++;
		}
	}


	buildContour_wavy(contourPanelSeparation, shapeHeight, shapeDepth, contours_vert);

	generateRandomDots(shapeWidth, shapeHeight, shapeDepth, dot_num_per_col, dot_num_per_col, dot_jitter_max_scale, dotContainer);

}

void generateRandomDots(double shapeWidth, double shapeHeight, double shapeDepth, int dotNumPerRow, int dotNumPerCol, double dotJitterMax_Scale, vector<Vector3d>& dotContainer) {

	dotContainer.clear();
	int dot_counter_top = 0, dot_counter_bottom = 0;
	int N_count = 2;

	double step_size_x = shapeWidth / (double)(dotNumPerRow - 1);
	double step_size_y = shapeHeight / (double)(dotNumPerCol - 1);

	double dot_jitter_x_max = step_size_x * dotJitterMax_Scale;
	double dot_jitter_y_max = step_size_y * dotJitterMax_Scale;

	double x, y, z;

	for (int i_y = 0; i_y < dotNumPerCol; i_y++) {
		for (int i_x = 0; i_x < dotNumPerRow; i_x++) {
			x = i_x * step_size_x - shapeWidth / 2 +
				(rand() % 17) / 16.0 * dot_jitter_x_max - dot_jitter_x_max / 2;

			y = i_y * step_size_y - shapeHeight / 2 +
				(rand() % 17) / 16.0 * dot_jitter_y_max - dot_jitter_y_max / 2;

			if (y < -shapeHeight / 2) {
				y = -shapeHeight / 2;
				z = 0;
				dot_counter_bottom++;

				if (dot_counter_bottom % N_count == 0) {
					dotContainer.push_back(Vector3d(x, y, z));
				}

			}
			else if (y > shapeHeight / 2) {
				y = shapeHeight / 2;
				z = 0;
				dot_counter_top++;

				if (dot_counter_top % N_count == 0) {
					dotContainer.push_back(Vector3d(x, y, z));
				}
			}
			else {

				z = getZ(shapeHeight, shapeDepth, y);
				dotContainer.push_back(Vector3d(x, y, z));
			}

		}
	}
}


void drawRandomDots(const vector<Vector3d>& dotContainer) {

	glPushMatrix();
	glTranslated(0, 0, .8);

	glColor3f(0.0f, 0.0f, 0.0f);

	for (int i = 0; i < int(dotContainer.size()); i++)
	{
		Vector3d dot_vector = dotContainer.at(i);
		double x_axis = dot_vector[0];
		double y_axis = dot_vector[1];
		double z_axis = dot_vector[2];

		glPushMatrix();
		glTranslated(x_axis, y_axis, z_axis);

		double dot_size = ratio_dotSize_distance * abs(display_distance + z_axis);
		glutSolidSphere(dot_size, dot_sides, dot_sides);
		glPopMatrix();
	}

	glPopMatrix();

}


/***** SOUND THINGS *****/



void initGLVariables()
{

}

void drawInfo()
{
	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_BLEND);

	GLText text;
	if (gameMode)
		text.init(SCREEN_WIDTH, SCREEN_HEIGHT, glWhite, GLUT_BITMAP_HELVETICA_18);
	else
		text.init(640, 480, glWhite, GLUT_BITMAP_HELVETICA_12);

	text.enterTextInputMode();
	glColor3fv(glRed);
	text.draw("#disparity depth: " + stringify<double>(depth_disp));
	text.draw("#visual_angle: " + stringify<double>(visual_angle));


	text.leaveTextInputMode();
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_BLEND);
}

void drawStimulus()
{
	if (stimulus_built) {
		drawRandomDotSurface(dist_toEye, my_verts, dot_container, my_contour_data);
	}
}


void drawGLScene()
{
	glDrawBuffer(GL_BACK_RIGHT);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0.0, 0.0, 0.0, 1.0);
	cam.setEye(eyeRight);
	//cam.setEye(eyeCyclopean);
	drawStimulus();
	drawInfo();

	// Draw right eye view
	glDrawBuffer(GL_BACK_LEFT);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0.0, 0.0, 0.0, 1.0);
	cam.setEye(eyeLeft);
	//cam.setEye(eyeCyclopean);
	drawStimulus();
	drawInfo();

	glutSwapBuffers();
	glutPostRedisplay();
}




// Funzione di callback per gestire pressioni dei tasti
void handleKeypress(unsigned char key, int x, int y)
{
	switch (key)
	{   //Quit program


	case 27:	//corrisponde al tasto ESC
	{
		if (resetScreen_betweenRuns)
			homeEverything(5000, 4500);

		exit(0);
	}
	break;

	// Enter key: press to make the final calibration
	case '+':
	{
		initStimulus();
	}
	break;

	case '1':
	{
		if (depth_disp > (5 + depth_inc)) {
			depth_disp = depth_disp - depth_inc;
			initStimulus();
		}

	}
	break;

	case '2':
	{
		depth_disp = depth_disp + depth_inc;
		initStimulus();
	}
	break;


	case '4':
	{

	}
	break;

	case '5':
	{

	}
	break;

	case '7':
	{

	}
	break;

	case '8':
	{

	}
	break;


	case '0':
		visual_angle = visual_angle - 0.2;
		initStimulus();
		break;

	case '.':
		visual_angle = visual_angle + 0.2;
		initStimulus();
		break;




	case '3':

		break;

	case '6':
		break;





	}
}


void handleResize(int w, int h)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
}


void initProjectionScreen(double _focalDist, const Affine3d& _transformation, bool synchronous)
{
	focalDistance = _focalDist;
	screen.setWidthHeight(SCREEN_WIDE_SIZE, SCREEN_WIDE_SIZE * SCREEN_HEIGHT / SCREEN_WIDTH);
	screen.setOffset(alignmentX, alignmentY);
	screen.setFocalDistance(_focalDist);
	screen.transform(_transformation);
	cam.init(screen);
	if (synchronous)
		moveScreenAbsolute(_focalDist, homeFocalDistance, 3500);
	else
		moveScreenAbsoluteAsynchronous(_focalDist, homeFocalDistance, 3500);
}


// Questa funzione si occupa di fare il refresh della schermata ed e' chiamata ogni TIMER_MS millisecond, tienila cosi'
void update(int value)
{
	glutPostRedisplay();
	glutTimerFunc(TIMER_MS, update, 0);
}

void idle()
{
	glutPostRedisplay();

}

void initRendering()
{
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	/* Set depth buffer clear value */
	glClearDepth(1.0);
	/* Enable depth test */
	glEnable(GL_DEPTH_TEST);
	/* Set depth function */
	glDepthFunc(GL_LEQUAL);
	// scommenta solo se vuoi attivare lo shading degli stimoli

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	// Tieni questa riga per evitare che con l'antialiasing attivo le linee siano piu' sottili di un pixel e quindi
	// ballerine (le vedi vibrare)
	glLineWidth(1.5);


}

void initVariables()
{
	eyeRight = Vector3d(interoculardistance / 2.0, 0, 0);
	eyeLeft = Vector3d(-interoculardistance / 2.0, 0, 0);
	eyeCyclopean = Vector3d(0, 0, 0);

	initProjectionScreen(display_distance);

}


// Inizializza gli stream, apre il file per poi scriverci
void initStreams()
{

}
// Porta tutti i motori nella posizione di home e azzera i contatori degli steps
void initMotors()
{
	if (resetScreen_betweenRuns)
		homeEverything(5000, 4500);
}

int main(int argc, char* argv[])
{
	mathcommon::randomizeStart();
	glutInit(&argc, argv);

	//glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STEREO);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_STEREO);

	if (gameMode == false)
	{
		glutInitWindowSize(640, 480);
		glutCreateWindow("EXP WEXLER");
		//glutFullScreen();
	}
	else
	{
		glutGameModeString("1024x768:32@85");
		glutEnterGameMode();
		glutFullScreen();
	}
	initMotors();
	initRendering();
	initGLVariables();
	initStreams();
	initVariables();
	initStimulus();
	glutDisplayFunc(drawGLScene);
	glutKeyboardFunc(handleKeypress);
	glutReshapeFunc(handleResize);
	glutIdleFunc(idle);
	glutTimerFunc(TIMER_MS, update, 0);
	glutSetCursor(GLUT_CURSOR_NONE);
	//boost::thread initVariablesThread(&initVariables);
	glutMainLoop();

	return 0;
}
