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

#include "SOIL.h"//Library for texture mapping

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
//TrialGenerator<double> trial;//if using staircase: 

ParametersLoader parameters_subj;
ParametersLoader parameters;

// paramters file directory and name

string parametersFileName_subj = experiment_directory + "trialPars/parameters_Subj.txt";

string parametersFileName = experiment_directory + "trialPars/parameters_Haptic_probe.txt";

// response file
ofstream responseFile;
string responseFile_headers = "subjName\tIOD\tblockN\ttrialN\tdisplayDistance\tvisualAngle\tshapeHeight\tshapeID\ttestingTexture\ttestingDisparity\tlefEyeView\trightEyeView\ttexnum\ttextNomralizer\ttestDepth\tprobeDepthInit\tprobeDepth\tdotNum\tdotVisAg\tdotColR\tjitterScale\tRT";

string subjectName;


/**********	TRIALS **************/
int sessionNum = 0;
int totalBlkNum = 1;
int blkNum = 1;
int trialNum = 0;

int trainNum_cap = 50;

double percentComplete = 0;
int repetition = 5;
int totalTrNum = 6 * repetition;

/********** STIMULUS SHAPE ***************/
// stimulus shape
double display_distance;
double visual_angle = 7.5; // stim diangonal size
double jitter_z = 0;
double display_distance_jittered = display_distance + jitter_z;

double stimulus_height = 70;//tan((DEG2RAD * visual_angle)/2) * 2 * (abs(display_distance));
double stimulus_width = 70;
double ratio_width_height = 1.06;
double ratio_width_height_long = 1.4;

double depth_test = 35;
double depth_inc = 2;
double depth_training_min = 2; // set in the subject file
int depth_training_range = 50; // set in the subject file

enum shapeTypes {Ridge, Gaussian, Cosine, CosineRidge };
shapeTypes current_shape = CosineRidge;
int shapeID = 0;
double gauss_sig_height_ratio = 0.16;

int nr_points = 201;

/********** BLOCKING PANELS ***************/
enum panelStates{no_aperture, black_aperture, red_aperture};
panelStates panel_state = black_aperture;

std::vector<Vector3d> vertContainer_Lcontour;
std::vector<Vector3d> vertContainer_Rcontour;
/********** RAMDOM DOTS ***************/
std::vector<Vector3d> dot_container;

// dot specs
int dot_number = 180; //num = 500, in the square which spans 7 degree in visual angle at 400mm display distance
// for ot size
double visAngle_dot = 0.20;//Bigger dot: 0.28; smaller dot: 0.14
double ratio_dotSize_distance = tan(DEG2RAD * visAngle_dot / 2.0);
float dot_color_R = 0.75; //Bigger dot: 0.64; smaller dot: 0.85
int dot_sides = 16; // Bigger dot: 24; smaller dot: 8

// for lattice dot array
double dot_jitter_max_scale = 0.80;
int dot_num_per_col = 9;
/********** STIMULUS VERTICES ***************/
// vectors storing vertices data
std::vector<GLfloat> vertices_vec;
std::vector<GLfloat> texcoors_vec;
std::vector<GLfloat> colors_vec;
std::vector<GLfloat> normals_vec;
std::vector<GLuint> indices_draw_triangle_vec;

std::vector<GLfloat> vertices_vec_new;
std::vector<GLfloat> colors_vec_new;
std::vector<GLuint> indices_draw_triangle_vec_new;

// texture map
int texnum = 50;
GLuint loaded_textures[51];
double normalizer_to_uv_base = 90;
double normalizer_to_uv = 90;
double u_offset = 0.05;
double v_offset = 0.05;

float amb_intensity = 0.5;
/********** PROBE ***************/
double probe_depth = 10; // the depth indicated by the probe
double probe_depth_init = 30; // starting value, should be randomized
double probe_height = 100; // the height is the length when the probe is a straight line
double probe_location_offset = -80;
double probe_jitter_z = 0;

double probe_y_array[201];
double probe_x_array[201];

/********** STATE VARIABLES ***************/
enum Stages { stimulus_preview, prep_trial, trial_fixate, trial_view, trial_respond, break_time, exp_completed };
Stages current_stage = stimulus_preview; // if just want to look at the stimuli, use the constant present stage

//state bool 
bool testing_texture_vs_disparity = true;
bool monocular_display = true;
bool left_eye_on = true;
bool right_eye_on = true;

bool training = true;
bool visibleInfo = true;

bool resetScreen_betweenRuns = false;

/********** TIME ***************/
// Timer variable, set for each trial 
Timer trial_timer;
double ElapsedTime;

double fixateTime = 800;
double viewTime = 800;
double progressBarTime = 1500;
double responseTime = 8000;

/*********** for DEBUGGING **********/


/*************************** FUNCTIONS ***********************************/
void initOptotrak();
void initMotors();
void initRendering();
void initVariables();
void initStreams();
void initBlock();
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

void buildVertices_textureMap(double shapeDepth);
std::vector<Vector3d> buildRandomDots_Lattice(double shapeDepth, int dotNumPerRow, int dotNumPerCol, double dotJitterMax_Scale);
void drawVertices_textureMap(int texNum, double dispDist, double shapeDepth);
std::vector<Vector3d> buildRandomDots(double shapeDepth);
void drawRandomDots(std::vector<Vector3d> dot_container, double dispDist, double shapeDepth);
void updateProbe(double probeHeight, double probeDepth);
void drawProbe(double dispDist, double probeDepth);
void drawProgressBar();
void drawTimeBar();
void drawBlockingPanels(double pandelSeparation);
void drawFixation(double dispDist);

int LoadGLTextures();
void setLightingPar();
void buildVertices_randomDot(double shapeDepth);
void drawVertices_randomDot(double dispDist, double shapeDepth);
void buildVertices_panels_wavy(double shapeDepth);
void buildVertices_panels_jigsaw(double shapeDepth);
void drawVertices_panels();
