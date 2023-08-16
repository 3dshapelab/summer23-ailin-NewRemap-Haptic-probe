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

/********** TEXTURE SURF ***************/
int nr_curve_map = 10001;

double lengthFactor_TM = 1.4; //for a depth with a curve length of l, the TM length is multiplied by this factor.
double del_l = 0.4;

int nr_points_width = 251; // nr of points in x direction
int nr_points_height_default = 201; // default
int nr_points_height = nr_points_height_default;
int total_ind = 0;

/********* TEXTURE *********/
// self-generated
float Tex_dot_radius = 2.7;
float Tex_dot_density = 0.019;
float Tex_dot_separation_ratio = 1.10;

float texture_col_max = 1.0;
float texture_col_min = 0.1;
int TexDot_Lat_nr = 4;
float TexDot_Lat_jitter = 0.4;

//blur edge
double drop_off_rate = 0.75;
double R_intersect_factor = 2 / (1 + drop_off_rate);


/********** LIGHT SHADING ***************/
float max_intensity = 1.0;
float light_amb = 0.3;
float light_dif = 0.5;
float lightDir_z = 0.5;
double light_depthMin = 26;
double light_depthMax = 40;


/********** PROBE ***************/
bool probe_drawSideView = true;
double probe_depth = 10; // the depth indicated by the probe
double probe_depth_init = 30; // starting value, should be randomized
double probe_height = 100; // the height is the length when the probe is a straight line
double probe_location_offset = -80;
double probe_jitter_z = 0;
int nr_points_probe = 101;
vector<Vec2> probe_verts;

/********** OPTIONS ***************/
bool resetScreen_betweenRuns = false;
enum surfaceTypes{randomdots_surface, polkadots_surface};
surfaceTypes current_surface = randomdots_surface; 

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

void generateRandomDots(double shapeWidth, double shapeHeight, double shapeDepth, int dotNumPerRow, int dotNumPerCol, double dotJitterMax_Scale, vector<Vector3d>& dotContainer);
void buildContour_wavy(double ContourWidth, double shapeHeight, double shapeDepth, ContourData& new_contours_vert);
void buildRandomDotSurface(double shapeWidth, double shapeHeight, double shapeDepth, double contourPanelSeparation, VerticesData& vertices_data, vector<Vector3d>& dotContainer, ContourData& contours_vert);
void drawRandomDots(const vector<Vector3d>& dotContainer);
void drawRandomDotSurface(double distShapeToEye, const VerticesData& vertices_data, const vector<Vector3d>& dotContainer, const ContourData& contours_vert);

void generateTexture(float TM_X, float TM_Y, float dotDensity, float dotRadius, float dotSeparationRatio_init, int nr_X_Lattice, float dotJitterScale_Lattice, TextureDotsData& outputTexDots);
void buildVertices_Texture(double shapeWidth, const CurvePtsData& dispYCurve, const CurvePtsData& textYCurve, double distShapeToEye, TextureDotsData& TexDotsOnText, VerticesData& vertices_data);
void buildContour_Texture(double ContourWidth, const CurvePtsData& dispYCurve, const CurvePtsData& textYCurve, float distShapeToEye, ContourData& new_contours_vert);
void buildTextureSurface(double shapeWidth, double shapeHeight, double dispDepth, double textDepth, double distShapeToEye, double contourPanelSeparation, VerticesData& vertices_data, ContourData& contours_vert);
void drawTextureSurface(double distShapeToEye, const VerticesData& vertices_data, const ContourData& contours_vert);
void drawContours(const ContourData& contours_vert);

double getZ(double shapeHeight, double shapeDepth, double vertexY);
double getTg(double shapeHeight, double shapeDepth, double Y);
double SolveForZ_projected(double theHeight, double newDepth, double l, double y0, double z0);
void scanCurve(double shapeHeight, double shapeDepth, CurveYLMap& output_curve_ylmap);
void projectCurve(const CurveYLMap& curve_map_proj, double distShapeToEye, const CurvePtsData& origin_curve, CurvePtsData& output_curve_proj);
Vector3d projectPoint(double shapeHeight, double newDepth, double distShapeToEye, Vector3d fromPoint);
float adjustDiffLight(double textDepth, float maxInt, float ambInt, double Depth_flat, double Depth_deep);

void updateProbe(double probeHeight, double probeDepth, vector<Vec2>& probeVerts);
void drawProbe(double dispDist, double probeDepth, const vector<Vec2>& probeVerts);

/*************************** FUNCTIONS: random dots surface ***********************************/

double getZ(double shapeHeight, double shapeDepth, double Y) {

	double Z;

	Z = shapeDepth * cos(M_PI * Y / shapeHeight);

	return (Z);
}

void updateProbe(double probeHeight, double probeDepth, vector<Vec2>& probeVerts) {

	probeVerts.clear();

	double step_size = probeHeight / (double)(nr_points_probe - 1);

	for (int j = 0; j < nr_points_probe; j++) {

		float y = -probeHeight / 2 + j * step_size;
		float x = getZ(probeHeight, probeDepth, y);
		probeVerts.push_back(Vec2{ x,y });

	}

}

void drawProbe(double dispDist, double probeDepth, const vector<Vec2>& probeVerts) {

	glDisable(GL_TEXTURE_2D);
	//glColor3f(1.0f, 0.0f, 0.0f);

	glPushMatrix();
	glLoadIdentity();
	glTranslated(0, probe_location_offset, dispDist);

	if (probe_drawSideView) {

		glBegin(GL_LINES);
		glLineWidth(1.f);
		glColor3f(0.5f, 0.0f, 0.0f);
		for (int i = 0; i < probeVerts.size() - 1; i++) {
			glVertex3d(probeVerts[i].x, probeVerts[i].y, 0);
			glVertex3d(probeVerts[i+1].x, probeVerts[i+1].y, 0);
		}

		glEnd();
	}

	glBegin(GL_LINES);
	glLineWidth(2.f);
	glColor3f(0.8f, 0.0f, 0.0f);
	glVertex3d(0, -40, 0);
	glVertex3d(0, 40, 0);
	glVertex3d(probeDepth, -40, 0);
	glVertex3d(probeDepth, 40, 0);
	glEnd();


	glPopMatrix();
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

/*************************** FUNCTIONS: texture surface ***********************************/
double getTg(double shapeHeight, double shapeDepth, double Y) {
	return (-shapeDepth * sin(M_PI * Y / shapeHeight) * M_PI / shapeHeight);
}

float adjustDiffLight(double textDepth, float maxInt, float ambInt, double Depth_flat, double Depth_deep) {
	float newDiff = ambInt;
	newDiff = newDiff + (textDepth - Depth_flat) / (Depth_deep - Depth_flat) * (maxInt - 2 * ambInt);
	return newDiff;
}

double NewtonSolver_fz(double z, double Depth, double zCoeff, double distShapeToEye) {
	double val = z / Depth - cos(zCoeff * (z - distShapeToEye));
	return val;
}

double NewtonSolver_dfz(double z, double Depth, double zCoeff, double distShapeToEye) {
	double val = 1 / Depth + sin(zCoeff * (z - distShapeToEye)) * zCoeff;
	return val;
}

double SolveForZ_projected(double theHeight, double newDepth, double distShapeToEye, double y0, double z0) {

	double z_new, z, f_z, df_z;
	double C = M_PI * y0 / (theHeight * (z0 - distShapeToEye));

	z = z0;

	for (int i = 0; i < 100; i++) {

		f_z = NewtonSolver_fz(z, newDepth, C, distShapeToEye);
		df_z = NewtonSolver_dfz(z, newDepth, C, distShapeToEye);

		if (abs(f_z) < 1e-10) {

			break;
		}
		else if (abs(df_z) < 1e-10) {

			break;
		}
		else {
			z_new = z - f_z / df_z;
			z = z_new;
		}
	}

	if (abs(z - z0) > 40)
		z = 0;

	return z;

}

Vector3d projectPoint(double shapeHeight, double newDepth, double distShapeToEye, Vector3d fromPoint) {

	Vector3d ToPoint = fromPoint;
	if (abs(abs(fromPoint.y()) - shapeHeight / 2) > 0.01) {
		double z = SolveForZ_projected(shapeHeight, newDepth, distShapeToEye, fromPoint.y(), fromPoint.z());
		double w = (z - distShapeToEye) / (fromPoint.z() - distShapeToEye);
		ToPoint = Vector3d(w * fromPoint.x(), w * fromPoint.y(), z);
	}

	return ToPoint;
}

double normalCDF(double value)
{
	double M_SQRT1_2 = sqrt(0.5);
	return 0.5 * erfc(-value * M_SQRT1_2);
}

double blurEdge(double dropoff, double pointDist, double intersectDist, double baseCol, double maxCol) {

	// given an Radius, I will choose the intersectDist to be 2/(1+dropoff)
	double addCol = maxCol - baseCol;
	return(((pointDist / intersectDist) < dropoff) ? baseCol : (baseCol + addCol * normalCDF(((pointDist / intersectDist - dropoff) / (1 - dropoff) * 6) - 3)));
}

void scanCurve(double shapeHeight, double shapeDepth, CurveYLMap& output_curve_ylmap) {

	// go through the cosine curve with input height and depth, equally sample nr_curve_map dots and mark the paired y_l values at these points
	output_curve_ylmap = {};

	double y, z, l, y_prev, z_prev;
	double step_size = (shapeHeight / (nr_curve_map - 1));

	output_curve_ylmap.curve_depth = shapeDepth;
	output_curve_ylmap.curve_height = shapeHeight;
	output_curve_ylmap.step_size = step_size;

	// the first point
	y = -shapeHeight / 2;
	z = 0;
	l = 0;
	y_prev = y, z_prev = z;
	output_curve_ylmap.y_vec.push_back(y);
	output_curve_ylmap.l_vec.push_back(l);


	for (int j = 1; j < nr_curve_map; j++) {
		y = -shapeHeight / 2 + j * step_size;
		z = getZ(shapeHeight, shapeDepth, y);
		l = l + sqrt(pow(y - y_prev, 2) + pow(z - z_prev, 2));

		output_curve_ylmap.y_vec.push_back(y);
		output_curve_ylmap.l_vec.push_back(l);

		y_prev = y; z_prev = z;
	}


}

void projectCurve(const CurveYLMap& curve_map_proj, double distShapeToEye, const CurvePtsData& origin_curve, CurvePtsData& output_curve_proj) {

	output_curve_proj = {};

	double newDepth = curve_map_proj.curve_depth;
	double height = curve_map_proj.curve_height;
	double step_y_ylmap = curve_map_proj.step_size;

	output_curve_proj.curve_height = height;
	output_curve_proj.curve_depth = newDepth;

	double y_p, z_p, l_p, tg_p;


	for (int jj = 0; jj < origin_curve.y_vec.size(); jj++) {

		double y_o = origin_curve.y_vec[jj];
		double z_o = origin_curve.z_vec[jj];
		z_p = SolveForZ_projected(height, newDepth, distShapeToEye, y_o, z_o);
		double w = (z_p - distShapeToEye) / (z_o - distShapeToEye);
		y_p = w * y_o;
		int i_c = (y_p - curve_map_proj.y_vec[0]) / step_y_ylmap;
		l_p = curve_map_proj.l_vec[i_c];

		output_curve_proj.y_vec.push_back(y_p);
		output_curve_proj.z_vec.push_back(z_p);
		output_curve_proj.l_vec.push_back(l_p);

	}

}

int buildCurve_byDelY(const CurveYLMap& input_curve_ylmap, CurvePtsData& output_curve) {

	output_curve = {};

	double depth = input_curve_ylmap.curve_depth;
	output_curve.curve_depth = depth;
	double height = input_curve_ylmap.curve_height;
	output_curve.curve_height = height;
	double y, l;

	double stpsz_J = height / (nr_points_height_default - 1);
	double stpsz_ycurve_precise = height / (nr_curve_map - 1);
	for (int j = 0; j < nr_points_height_default; j++) {
		int k = stpsz_J * j / stpsz_ycurve_precise;
		y = input_curve_ylmap.y_vec[k];
		output_curve.y_vec.push_back(y);
		output_curve.z_vec.push_back(getZ(height, depth, y));
		output_curve.l_vec.push_back(input_curve_ylmap.l_vec[k]);
	}

	return output_curve.y_vec.size();
}

void generateTexture(float TM_X, float TM_Y, float dotDensity, float dotRadius, float dotSeparationRatio_init, int nr_X_Lattice, float dotJitterScale_Lattice, TextureDotsData& outputTexDots) {

	// the hybrid here refers to employing two methods to generate dots: 
	// method 1: lattice with some jitter
	// method 2: random placement with minimum distance constraint
	// 
	// INPUT:
	// TM_X & TM_Y control the dimension of texture map
	// dotDensity controls how dense the dots are
	// dotSeparationRatio_init controls the minimum distance between two dots, but it is an initial value, it may decrease if it constrains placing all dots in desingiated area
	// nr_X_Lattice controls the number of dots on each row of the lattice
	// dotJitterScale_Lattice controls the scale of jitter of dots on lattice
	// OUTPUT:
	// outputTexDots are TextureDotsData, the most important component is a vector of dot centers (coordinate origin is at the center of the TM) 
	// the dot center x takes value from -TM_X/2 to TM_X/2, the dot center y takes value from -TM_Y/2 to TM_Y/2

	outputTexDots = {};

	outputTexDots.TexMapSize = Vec2{ TM_X, TM_Y };
	outputTexDots.Radius = dotRadius;
	outputTexDots.margin_y = R_intersect_factor * dotRadius;

	std::uniform_real_distribution<float> dist(0.f, 1.f); // require <random>


	vector<Vec2> dc_vec;

	int num_dot = TM_X * TM_Y * dotDensity;

	// prep lattice
	int nr_X = nr_X_Lattice;
	int nr_Y = floor(TM_Y * nr_X / TM_X);
	if (nr_Y % 2 == 0) {
		nr_Y++;
	}
	float lat_step_x = TM_X / (float)nr_X;
	float lat_step_y = TM_Y / (float)nr_Y;



	// step 1: generate lattice dots
	for (int j = 0; j < nr_Y; j++) {
		for (int i = 0; i < nr_X; i++) {
			float rx = dist(rng);
			float ry = dist(rng);
			float cx_lat = ((rx - 0.5) * dotJitterScale_Lattice + i + 0.5) * lat_step_x;
			float cy_lat = ((ry - 0.5) * dotJitterScale_Lattice + j + 0.5) * lat_step_y;

			dc_vec.push_back(Vec2{ cx_lat, cy_lat });
		}
	}


	// step 2: random placement with minimum distance constraint
	float dotSeparationRatio_adjust = dotSeparationRatio_init;
	float dot_separation = dotSeparationRatio_adjust * 2 * (dotRadius * R_intersect_factor);

	int reset_count = 0;
	int dotplacement_runs = 0;

	for (int i_dot = nr_X * nr_Y; i_dot < num_dot; i_dot++) {

		dotplacement_runs++;

		// step 2.1: check whether needs to reset or exit
		// too many unsuccessful placements leads to reset, too many resets leads to exit()
		if (dotplacement_runs > 10000) {

			reset_count++;

			if (reset_count < 10000) {

				//too many unsuccessful placements leads to reset
				dotplacement_runs = 0;
				dc_vec.clear();

				for (int j = 0; j < nr_Y; j++) {
					for (int i = 0; i < nr_X; i++) {

						float rx = dist(rng);
						float ry = dist(rng);
						float cx_lat = ((rx - 0.5) * dotJitterScale_Lattice + i + 0.5) * lat_step_x;
						float cy_lat = ((ry - 0.5) * dotJitterScale_Lattice + j + 0.5) * lat_step_y;

						dc_vec.push_back(Vec2{ cx_lat, cy_lat });
					}
				}

				i_dot = nr_X * nr_Y;

				int adjust_count = floor(reset_count / (float)50); // make one adjustment after repeating fifty attempts
				dotSeparationRatio_adjust = dotSeparationRatio_adjust - 0.000005 * (float)adjust_count * (float)num_dot;
				dot_separation = dotSeparationRatio_adjust * 2 * (dotRadius * R_intersect_factor);
			}
			else {
				// too many resets leads to exit()
				cout << "text depth: " << depth_text << endl;
				cout << "disp depth: " << depth_disp << endl;
				cerr << "can't generate texture" << endl;
				exit(0);
			}

		}


		// step 2.2: placing one dot at a time
		// pick random xy values for the circle center
		float cx = dist(rng); // 0-1
		float cy = dist(rng); // 0-1
		cx *= TM_X;
		cy *= TM_Y;

		// checking whether the current circle intersects withe the previous circles already pushed
		bool intersect = false;
		for (int k = 0; k < dc_vec.size(); k++) {
			Vec2 prev_c = dc_vec[k];
			if ((cx - prev_c.x) * (cx - prev_c.x) + (cy - prev_c.y) * (cy - prev_c.y) < dot_separation * dot_separation) {
				intersect = true;
				break;
			}
		}
		// if intersect, then break and pick a new circle center
		if (intersect) {
			i_dot--;
			continue;
		}

		// if not intersect, add this circle to the circles vector
		dc_vec.push_back(Vec2{ cx, cy });
	}


	// step 3: sort the dots by their y. Get the correponding indices
	int n_dots = dc_vec.size();

	vector<float> dc_y_vec;
	for (int k = 0; k < n_dots; k++) {
		dc_y_vec.push_back(dc_vec[k].y);
	}
	vector<int> sortedDC_ind_vec(n_dots);
	std::iota(sortedDC_ind_vec.begin(), sortedDC_ind_vec.end(), 0); //Initializing
	std::sort(sortedDC_ind_vec.begin(), sortedDC_ind_vec.end(), [&](int i, int j) {return dc_y_vec[i] < dc_y_vec[j]; });



	// step 4: fill in outputTexDots dot center vectors, in the order of from low y to high y
	float x_offset = TM_X / 2.;
	float y_offset = TM_Y / 2.;
	for (int k = 0; k < n_dots; k++) {
		Vec2 dc_temp = dc_vec[sortedDC_ind_vec[k]];
		outputTexDots.dot_center_vec.push_back(Vec2{ dc_temp.x - x_offset, dc_temp.y - y_offset });
	}

}

void buildContour_Texture(double ContourWidth, const CurvePtsData& dispYCurve, const CurvePtsData& textYCurve, float distShapeToEye, ContourData& new_contours_vert) {

	new_contours_vert = {};

	for (int i_v = 0; i_v < dispYCurve.y_vec.size(); i_v++) {

		float x_t_L = -ContourWidth / 2;
		float x_t_R = ContourWidth / 2;

		float z_t = textYCurve.z_vec[i_v];

		float y_d = dispYCurve.y_vec[i_v];
		float z_d = dispYCurve.z_vec[i_v];

		float w = (distShapeToEye - z_d) / (distShapeToEye - z_t);
		//float w = (distShapeToEye - z_d) / (distShapeToEye - (z_t + z_d) / 2.0);
		float x_v_L = w * x_t_L;
		float x_v_R = w * x_t_R;

		new_contours_vert.vert_Lcontour.push_back(Vector3f(x_v_L, y_d, z_d));
		new_contours_vert.vert_Rcontour.push_back(Vector3f(x_v_R, y_d, z_d));
	}
}

void buildVertices_Texture(double shapeWidth, const CurvePtsData& dispYCurve, const CurvePtsData& textYCurve, double distShapeToEye, TextureDotsData& TexDotsOnText, VerticesData& vertices_data) {

	vertices_data = {};

	GLuint i_ind = 0;
	int nr_J = dispYCurve.y_vec.size();
	double stpsz_I = shapeWidth / (nr_points_width - 1);
	double height = textYCurve.curve_height;
	double depth_text = textYCurve.curve_depth;
	float x_offset = TexDotsOnText.TexMapSize.x / 2;

	// for texture colors
	float vertex_col = 1.0f;
	int nr_dots = TexDotsOnText.dot_center_vec.size();
	float R = TexDotsOnText.Radius;

	float L_start = -textYCurve.l_vec.back() / 2.;
	float l_margin = TexDotsOnText.margin_y;

	std::uniform_real_distribution<float> dist(0.f, 1.f); // require <random>
	float L_offset = dist(rng) - 0.5f;
	//L_offset = L_offset * (TexDotsOnText.TexMapSize.y - textYCurve.l_vec.back() - 2 * TexDotsOnText.margin_y);
	L_offset = L_offset * 2 * R;
	L_start = L_start + L_offset;

	int TexDot_Ind_L = 0;
	int TexDot_Ind_H = 0;
	vector<Vec2> nearTexDots_dc_vec;

	for (int jj = 0; jj < nr_J; jj++) {

		float y_d = dispYCurve.y_vec[jj];
		float z_d = dispYCurve.z_vec[jj];

		float y_t = textYCurve.y_vec[jj];
		float z_t = textYCurve.z_vec[jj];

		float tg_t = getTg(height, depth_text, y_t);

		float w = (distShapeToEye - z_d) / (distShapeToEye - z_t);
		float x_d;

		double TM_y_t = textYCurve.l_vec[jj] + L_start;

		// find the dots that are near TM_y_t
		while (((TexDotsOnText.dot_center_vec[TexDot_Ind_L].y) < ((float)TM_y_t - l_margin)) && TexDot_Ind_L < (nr_dots - 1)) {
			TexDot_Ind_L++;
		}


		TexDot_Ind_H = TexDot_Ind_L;
		while (((TexDotsOnText.dot_center_vec[TexDot_Ind_H].y) < ((float)TM_y_t + l_margin)) && TexDot_Ind_H < (nr_dots - 1)) {
			TexDot_Ind_H++;
		}

		nearTexDots_dc_vec.clear();
		for (int k = TexDot_Ind_L; k < TexDot_Ind_H + 1; k++) {
			nearTexDots_dc_vec.push_back(TexDotsOnText.dot_center_vec[k]);
		}
		int nr_dots_near = nearTexDots_dc_vec.size();


		for (int ii = 0; ii < nr_points_width; ii++) {

			double pt_x = stpsz_I * ii - x_offset;
			double pt_y = TM_y_t;
			vertex_col = 1.0f;



			for (int k = 0; k < nr_dots_near; k++) {

				double pt_val = sqrt((pt_x - nearTexDots_dc_vec[k].x) * (pt_x - nearTexDots_dc_vec[k].x) +
					(pt_y - nearTexDots_dc_vec[k].y) * (pt_y - nearTexDots_dc_vec[k].y));


				if (pt_val < R_intersect_factor * R) {
					double vertex_col_tentative = blurEdge(drop_off_rate, pt_val, R_intersect_factor * R, texture_col_min, texture_col_max);
					vertex_col = float(vertex_col_tentative);
					break;
				}


				/*
								// not bluring the dot
								if (pt_val < R) {
									vertex_col = 0.0;
									break;
								}
				*/
			}

			x_d = w * pt_x;


			vertices_data.vertices_vec.push_back(x_d);
			vertices_data.vertices_vec.push_back(y_d);
			vertices_data.vertices_vec.push_back(z_d);

			vertices_data.light_normals_vec.push_back(0);
			vertices_data.light_normals_vec.push_back(-tg_t);
			vertices_data.light_normals_vec.push_back(1);


			vertices_data.colors_vec.push_back(vertex_col);
			vertices_data.colors_vec.push_back(0);
			vertices_data.colors_vec.push_back(0);

			if (ii < nr_points_width - 1 && jj < nr_J - 1) {

				// using vector
				vertices_data.indices_draw_triangle_vec.push_back(i_ind);
				vertices_data.indices_draw_triangle_vec.push_back(i_ind + 1);
				vertices_data.indices_draw_triangle_vec.push_back(i_ind + nr_points_width);

				vertices_data.indices_draw_triangle_vec.push_back(i_ind + nr_points_width);
				vertices_data.indices_draw_triangle_vec.push_back(i_ind + 1);
				vertices_data.indices_draw_triangle_vec.push_back(i_ind + nr_points_width + 1);

			}

			i_ind++;
		}
	}

}

void buildTextureSurface(double shapeWidth, double shapeHeight, double dispDepth, double textDepth, double distShapeToEye, double contourPanelSeparation, VerticesData& vertices_data, ContourData& contours_vert)
{

	// part 1: generate TexDots and project to Disp Surface
	CurveYLMap ylMap_Text;
	scanCurve(shapeHeight, textDepth, ylMap_Text);
	float l_text = ylMap_Text.l_vec.back();


	// part 3: generate surface vertices
	CurvePtsData y_curve_data_text_m, y_curve_data_disp_m;
	nr_points_height = buildCurve_byDelY(ylMap_Text, y_curve_data_text_m);

	if (abs(dispDepth - textDepth) < 0.1) {
		y_curve_data_disp_m = y_curve_data_text_m;
	}
	else {
		CurveYLMap ylMap_Disp;
		scanCurve(shapeHeight, dispDepth, ylMap_Disp);
		projectCurve(ylMap_Disp, distShapeToEye, y_curve_data_text_m, y_curve_data_disp_m);
	}


	//Tex_dot_radius = Tex_dot_radius_base + (float)(rand() % 11) * 0.012;
	TextureDotsData Tex_Dots_text;
	generateTexture(shapeWidth, lengthFactor_TM * l_text, Tex_dot_density, Tex_dot_radius, Tex_dot_separation_ratio, TexDot_Lat_nr, TexDot_Lat_jitter, Tex_Dots_text);

	buildVertices_Texture(shapeWidth, y_curve_data_disp_m, y_curve_data_text_m, distShapeToEye, Tex_Dots_text, vertices_data);

	buildContour_Texture(contourPanelSeparation, y_curve_data_disp_m, y_curve_data_text_m, distShapeToEye, contours_vert);

}

void drawTextureSurface(double distShapeToEye, const VerticesData& vertices_data, const ContourData& contours_vert) {

	//setting the light
	glShadeModel(GL_SMOOTH); // enable Smooth Shading
	glEnable(GL_LIGHTING); // enable lighting
	glEnable(GL_LIGHT1);
	glEnable(GL_NORMALIZE); //so we don't need to normalize our normal for surfaces	

	// Light source parameters
	GLfloat LightAmbient[] = { light_amb, 0.0f, 0.0f, 1.0f }; // non-directional & overall light (r,g,b,alpha): dark part
	GLfloat LightDiffuse[] = { light_dif, 0.0f, 0.0f, 1.0f }; // light created by the light source (directional light; r,g,b,alpha): bright part
	GLfloat LightPosition[] = { 0.0f, 1.f, lightDir_z, 0.0f }; // Light Position (x, y, z, 1.0f); if w==0, directional; if w==1, positional lights. Attenuation can be applied only to the positional light 

	//setting the light
	glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient); //setup the ambient light
	glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse); //setup the diffuse light
	glLightfv(GL_LIGHT1, GL_POSITION, LightPosition); //position the light

	glPushMatrix();
	glLoadIdentity();
	glTranslated(0, 0, -distShapeToEye);

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);

	// activate and specify pointer to vertex array
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	//using vector
	glVertexPointer(3, GL_FLOAT, 0, &vertices_data.vertices_vec[0]);
	glNormalPointer(GL_FLOAT, 0, &vertices_data.light_normals_vec[0]); //
	glColorPointer(3, GL_FLOAT, 0, &vertices_data.colors_vec[0]);
	glDrawElements(GL_TRIANGLES, vertices_data.indices_draw_triangle_vec.size(), GL_UNSIGNED_INT, &vertices_data.indices_draw_triangle_vec[0]);

	// deactivate vertex arrays after drawing
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);

	glDisable(GL_LIGHTING);

	drawContours(contours_vert);
	glPopMatrix();

}

/*************************** FUNCTIONS: general trial functions ***********************************/
void initStimulus() {

	stimulus_built = false;

	switch (current_surface)
	{
	case randomdots_surface:
		{
		dist_toEye = -(display_distance - depth_disp);
		stimulus_height = tan((DEG2RAD * visual_angle) / 2) * 2 * dist_toEye;
		stimulus_width = ratio_width_height_long * stimulus_height;
		stimulus_visiblewidth = ratio_visiblewidth_height * stimulus_height;
		buildRandomDotSurface(stimulus_width, stimulus_height, depth_disp, stimulus_visiblewidth, my_verts, dot_container, my_contour_data);
		}
		break;

	case polkadots_surface:
	{
		dist_toEye = -(display_distance - depth_disp);
		stimulus_height = tan((DEG2RAD * visual_angle) / 2) * 2 * dist_toEye;
		stimulus_width = ratio_width_height * stimulus_height;
		stimulus_visiblewidth = ratio_visiblewidth_height * stimulus_height;
		buildTextureSurface(stimulus_width, stimulus_height, depth_disp, depth_text, dist_toEye, stimulus_visiblewidth, my_verts, my_contour_data);
		light_dif = adjustDiffLight(depth_text, max_intensity, light_amb, light_depthMin, light_depthMax);
	}
	break;
	}
	
	updateProbe(stimulus_height, probe_depth, probe_verts);

	stimulus_built = true;
}


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
	text.draw("#texture depth: " + stringify<double>(depth_text));
	//text.draw("#visual_angle: " + stringify<double>(visual_angle));

	text.leaveTextInputMode();
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_BLEND);
}

void drawStimulus()
{
	if (stimulus_built) {

		switch (current_surface)
		{
		case randomdots_surface:
			drawRandomDotSurface(dist_toEye, my_verts, dot_container, my_contour_data);
		break;

		case polkadots_surface:		
			drawTextureSurface(dist_toEye, my_verts, my_contour_data);		
		break;
		}	


		drawProbe(display_distance + probe_jitter_z, probe_depth, probe_verts);
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
		current_surface = (surfaceTypes)((current_surface + 1) % 2 );
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
		if (depth_text > (5 + depth_inc)) {
			depth_text = depth_text - depth_inc;
			initStimulus();
		}
	}
	break;

	case '5':
	{
		depth_text = depth_text + depth_inc;
		initStimulus();
	}
	break;

	case '7':
	{
		probe_depth = probe_depth - 1;
		updateProbe(stimulus_height, probe_depth, probe_verts);
	}
	break;

	case '8':
	{
		probe_depth = probe_depth + 1;
		updateProbe(stimulus_height, probe_depth, probe_verts);
	}
	break;


	case '0':

		break;

	case '.':

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
