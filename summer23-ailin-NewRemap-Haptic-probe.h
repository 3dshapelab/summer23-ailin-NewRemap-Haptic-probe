#include "targetver.h"

#include <stdio.h>
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
#include <algorithm>
#include <queue>
#include <string>


/**** BOOOST MULTITHREADED LIBRARY *********/
#include <boost/thread/thread.hpp>
#include <boost/asio.hpp>	//include asio in order to avoid the "winsock already declared problem"


#ifdef _WIN32
#include <windows.h>
#include "CShaderProgram.h" // shader program
#include "glext.h" // opengl extensions for multisampling
#include <gl\gl.h>            // Header File For The OpenGL32 Library
#include <gl\glu.h>            // Header File For The GLu32 Library
#include "glut.h"            // Header File For The GLu32 Library
#include <MMSystem.h>
#endif

/************ INCLUDE CNCSVISION LIBRARY HEADERS ****************/
#include "Mathcommon.h"
#include "GLUtils.h"
#include "VRCamera.h"
#include "CoordinatesExtractor.h"
#include "GLText.h"

#include "ParametersLoader.h"
#include "Util.h"
#include "VRCamera.h"
#include "BalanceFactor.h"
#include "ParStaircase.h"
#include "Staircase.h"
#include "TrialGenerator.h"
#include "BrownPhidgets.h"
#include <direct.h>
#include "Optotrak2.h"
#include "Marker.h"
#include "BrownMotorFunctions.h"

#include <random>
#include <functional> // to help with thread safe randomization thread


/********* NAMESPACE DIRECTIVES ************************/
using namespace std;
using namespace mathcommon;
using namespace Eigen;
using namespace util;
using namespace BrownMotorFunctions;
using namespace BrownPhidgets;


/*************** Variable Declarations ********************/
static const bool gameMode = true;
const float DEG2RAD = M_PI / 180;

/********* VARIABLES OBJECTS  **********************/
VRCamera cam;
Optotrak2 optotrak;
Screen screen;
CoordinatesExtractor headEyeCoords;

/***** CALIBRATION FILE *****/
#include "Calibration_017B.h"
static const Vector3d center(0, 0, focalDistance);
double mirrorAlignment = 0.0, screenAlignmentY = 0.0, screenAlignmentZ = 0.0;

/********** EYES AND MARKERS **********************/
Vector3d eyeLeft, eyeRight, eyeMiddle;
vector <Marker> markers;
double interoculardistance = 60.0;
int screen1 = 19, screen2 = 20, screen3 = 21;
int mirror1 = 6, mirror2 = 22;


//////////////////////////////// usually no change is needed until this point //////////////////////


/*************************** TRIAL, INPUT AND OUTPUT ****************************/
// experiment directory
string experiment_directory = "C:/Users/labdomin/Documents/data/ailin/summer23-ailin-NewRemap-Haptic/";
BalanceFactor<double> trial; //if using costant stimuli

ParametersLoader parameters_subj;
ParametersLoader parameters;

// paramters file directory and name

string parametersFileName_subj = experiment_directory + "ParametersFiles/Haptic_Subj.txt";
string parametersFileName = experiment_directory + "ParametersFiles/parameters_Haptic_probe.txt";

// response file
ofstream responseFile;
string responseFile_headers = "subjName\ttrainCue\tIOD\tblockN\ttrialN\tdisplayDistance\tvisualAngle\ttestingTexture\ttestingDisparity\tlefEyeView\trightEyeView\tshapeHeight\tshapeWidth\tshapeDepth\tprobeDepthInit\tprobeDepth\tadjustUpNum\tadjustDnNum\tRT\tSurf_color\tnum_rdmDots\tvisAg_rdmDots\tnum_texDots\tradius_texDots";

string subjectName;


/**********	TRIALS **************/
int training_cue;
int sessionNum = 0;
bool sessionOrder_texture_first = false;
int totalBlkNum = 2;
int blkNum = 1;
int trialNum = 0;
int repetition = 2; // each block
int totalTrNum = totalBlkNum * repetition * 5; // two blocks * repetition * stimuli num
double percentComplete = 0;

int trainNum_cap = 8;
int adjustUpNum = 0, adjustDownNum = 0;

/********** STIMULUS SHAPE ***************/
// stimulus shape
double display_distance;
double visualTarget_X = 21;
double visual_angle = 7.4; // stim diangonal size
double jitter_z = 0;
double display_distance_jittered = display_distance + jitter_z;
double dist_toEye;

//height and width
double stimulus_height = 70; //tan((DEG2RAD * visual_angle)/2) * 2 * (abs(display_distance));
double stimulus_width = 70; //ratio_bgwidth_height * stimulus_height;
double stimulus_visiblewidth = 70; //ratio_visiblewidth_height * stimulus_height;
double ratio_width_height = 1.36;//1.3;
double ratio_visiblewidth_height = 1.15;//1.1;

double depth_test = 32;
double depth_inc = 2;
double depth_training_min = 20; // set in the subject file
int depth_training_range = 25; // set in the subject file

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
float dot_number;

/********** RAMDOM DOTS SURF ***************/
// dot distribution
int dot_num_per_col = 18;
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

/********** STATE VARIABLES ***************/
enum Stages { stimulus_preview, prep_trial, trial_fixate, trial_view, trial_respond, trial_confirm, trial_error, break_time, exp_completed };
Stages current_stage = stimulus_preview; // if just want to look at the stimuli, use the constant present stage

//state bool 
bool testing_texture_vs_disparity = true;
bool monocular_display = true;
bool left_eye_on = true;
bool right_eye_on = true;

bool training = true;
bool visibleInfo = true;
bool stimulus_built;


/********** TIME ***************/
// Timer variable, set for each trial 
Timer trial_timer;
double ElapsedTime;

double fixateTime = 800;
double viewTime = 800;
double progressBarTime = 1500;
double responseTime = 8000;

/********** OPTIONS ***************/
bool resetScreen_betweenRuns = false;

/*********** for DEBUGGING **********/


/*************************** FUNCTIONS ***********************************/
void initOptotrak();
void initMotors();
void initRendering();
void initVariables();
void initStreams();
void initBlock();
void setViewingVar();
void handleResize(int w, int h);
void initProjectionScreen(double _focalDist, const Affine3d& _transformation = Affine3d::Identity(), bool synchronous = true);
void updateTheMarkers();
void online_apparatus_alignment();
void cleanup();
void beepOk(int tone);

void initTrial();
void onlineTrial();
void advanceTrial();

void initStimulus(double stimulusDepth);
void drawStimulus();
void drawInfo();

void generateRandomDots(double shapeWidth, double shapeHeight, double shapeDepth, int dotNumPerRow, int dotNumPerCol, double dotJitterMax_Scale, vector<Vector3d>& dotContainer);
void buildContour_wavy(double ContourWidth, double shapeHeight, double shapeDepth, ContourData& new_contours_vert);
void buildRandomDotSurface(double shapeWidth, double shapeHeight, double shapeDepth, double contourPanelSeparation, VerticesData& vertices_data, vector<Vector3d>& dotContainer, ContourData& contours_vert);
void drawRandomDots(const vector<Vector3d>& dotContainer);
void drawRandomDotSurface(double distShapeToEye, const VerticesData& vertices_data, const vector<Vector3d>& dotContainer, const ContourData& contours_vert);

bool generateTexture(float TM_X, float TM_Y, float dotDensity, float dotRadius, float dotSeparationRatio_init, int nr_X_Lattice, float dotJitterScale_Lattice, TextureDotsData& outputTexDots);
void buildVertices_Texture(double shapeWidth, const CurvePtsData& dispYCurve, const CurvePtsData& textYCurve, double distShapeToEye, TextureDotsData& TexDotsOnText, VerticesData& vertices_data);
void buildContour_Texture(double ContourWidth, const CurvePtsData& dispYCurve, const CurvePtsData& textYCurve, float distShapeToEye, ContourData& new_contours_vert);
bool buildTextureSurface(double shapeWidth, double shapeHeight, double dispDepth, double textDepth, double distShapeToEye, double contourPanelSeparation, VerticesData& vertices_data, ContourData& contours_vert);
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

