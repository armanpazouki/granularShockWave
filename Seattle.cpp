// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Michal Kwarta
// =============================================================================
//
//
// =============================================================================

#include <iostream>
#include <vector>
#include <valarray>
#include <string>
#include <sstream>
#include <cmath>

#include "core/ChFileutils.h"
#include "core/ChStream.h"

#include "chrono_utils/ChUtilsCreators.h"
#include "chrono_utils/ChUtilsGenerators.h"
#include "chrono_utils/ChUtilsInputOutput.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::flush;
using std::endl;

// =============================================================
// Above headers and namespaces

//==============================================================
// 01. Velocity or force excitation?
// 02. DEM-P or DEM-C?
// 03. OpenGL?
// 04. Output directories and files 
// 05. Parameters describing bodies + gravity
// 06. Simulation time
// 07. Solver's parameters
// 08. Excitation's parameters
// 09. Output data's parameters
//==============================================================
// Functions:
//   I. AddWall
//  II. AddBall
// III. AddBalls
//==============================================================
// main(){...}

// =============================================================
// 01. Velocity or force excitation ?
// Decide what type of excitation you want to apply to the wall
//#define VELOCITY_EXCITATION

// =============================================================
// 02. DEM-P or DEM-C?
// Decide if you want to calculate with DEM-P or DEM-C
#define USE_DEM

// =============================================================
// 03. OpenGL?
// Decide if you want to see your simulation in OpenGL
#define CHRONO_PARALLEL_HAS_OPENGL

#ifdef CHRONO_PARALLEL_HAS_OPENGL
	#include "chrono_opengl/ChOpenGLWindow.h"
#endif

//======================================================
// 04. Output directories and files 
// Names of output directories and files
#ifdef USE_DEM
const std::string out_dir = "../Seattle_DEM";
#else
const std::string out_dir = "../Seattle_DVI";
#endif

const std::string pov_dir = out_dir + "/POVRAY";
const std::string forceSensor_file = out_dir + "/forceSensor.dat";

// -----------------------------------------------------------------------------
// Generate postprocessing output with current system state.
// -----------------------------------------------------------------------------
void OutputData(ChSystemParallel* sys, int out_frame, double time) {
  char filename[100];
  sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame);
  utils::WriteShapesPovray(sys, filename, false);
  cout << "time = " << time << flush << endl;
}

// =============================================================
// 05. Parameters describing bodies + gravity

double gravity = 9.81;

// spheres' dimentions & material
double ballRad = 0.5; // micrometer
double ballDiam = 2 * ballRad;

float mu = 0.2f;
float COR = 0.9f;
float ballDensity = 0.1f;
float ballMass = 1.0f;
float Y = 8e7; // to watch out for crazy vibrations - Y is 1000 times smaller
float nu = 0.3f;
ChVector<> inertia = (2.0 / 5.0) * ballMass * ballRad * ballRad * ChVector<>(1, 1, 1);

// container's dimentions & material
// number of spheres layers along x, y, z - axis
double xLay = 10;
double yLay = 21; // only even numbers - to hold symetry 
double zLay = 9;

double xDim = ballDiam * xLay;
double yDim = ballDiam + (yLay - 1.0) * sqrt(3.0) / 2.0 * ballDiam;
double zDim = ballDiam + (zLay - 1.0) * 2.0 * sqrt(2.0) / 3.0 * ballDiam;
double thickness = ballDiam;

float mu_ext = 0.2f; 
float COR_ext = 0.9f;
float density_ext = 0.1f;
float mass_ext = ballMass;
float Y_ext = 8e7; // to watch out for crazy vibrations - Y is 1000 times smaller
float nu_ext = 0.3f;

// =============================================================
// 06. Simulation time
// Parameters that describe the simulation
double simulation_time = 0.1;
// timeStartPushing (see Excitation parameters) is the time during 
// spheres settle themselves after it settled = true and wall is pushed
bool settled = false; 

// =============================================================
// 07. Solver's parameters

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 8;

// Perform dynamic tuning of number of threads?
bool thread_tuning = false;

// Solver settings
#ifdef USE_DEM
double time_step = 1.0e-4;
double tolerance = 0.1;
int max_iteration_bilateral = 1000;
#else
double time_step = 1e-3;
double tolerance = 0.1;
int max_iteration_normal = 0;
int max_iteration_sliding = 1000;
int max_iteration_spinning = 0;
int max_iteration_bilateral = 100;
double contact_recovery_speed = 10e30;
#endif
bool clamp_bilaterals = false;
double bilateral_clamp_speed = 0.1;

// =============================================================
// 08. Excitation's parameters
double timeStartPushing = 0.002; // start time for pushing the wall

#ifdef VELOCITY_EXCITATION
double timePushing = 0.0025; // period of time that the wall is pushed
double velocity = 500.0; // velocity of wall
#else
bool pushed = false; // to check if wall_5 came back a little bit to its initial position
double timePushing = time_step; // period of time that the wall is pushed
	#ifdef USE_DEM
	double force = 1e8;
	#else
	double force = 1e6;
	#endif
#endif

// =============================================================
// 09. Output data's parameters

bool write_povray_data = true;

double data_out_step = 1e-3;       // time interval between data outputs
double visual_out_step = 1e-2;     // time interval between PovRay outputs

int data_out_frame = 0;
int visual_out_frame = 0;

// =============================================================
// I. AddWall
// Function that adds the wall to the system
void AddWall(ChSystemParallel* sys, ChSharedPtr<ChBody> wall,
	int wallId, ChVector<> pos, ChVector<> dim,
	bool visualization){

	wall->SetIdentifier(wallId);
	//double wallMass = dim.x * dim.y * dim.z * wallDensity;
	wall->SetMass(mass_ext);
	wall->SetPos(pos);
	wall->SetBodyFixed(true);
	wall->SetCollide(true);

#ifdef USE_DEM
	ChSharedPtr<ChMaterialSurfaceDEM> mat_ext;
	mat_ext = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
	mat_ext->SetYoungModulus(Y_ext);
	mat_ext->SetPoissonRatio(nu_ext);
	mat_ext->SetRestitution(COR_ext);
	mat_ext->SetFriction(mu_ext);
#else
	ChSharedPtr<ChMaterialSurface> mat_ext;
	mat_ext = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	mat_ext->SetRestitution(COR_ext);
	mat_ext->SetFriction(mu_ext);
#endif

	wall->SetMaterialSurface(mat_ext);

	wall->GetCollisionModel()->ClearModel();

	utils::AddBoxGeometry(wall.get_ptr(), dim, ChVector<double>( 0, 0, 0), ChQuaternion<double>(1, 0, 0, 0), visualization);
	wall->GetCollisionModel()->SetFamily(1);
	wall->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

	wall->GetCollisionModel()->BuildModel();
	sys->AddBody(wall);
}

// =============================================================
// II. AddBall
// Function that adds the ball to the system

void AddBall(ChSystemParallel* sys, double x, double y, double z,	int ballId ){

#ifdef USE_DEM
	ChSharedPtr<ChMaterialSurfaceDEM> ballMat;
	ballMat = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
	ballMat->SetYoungModulus(Y);
	ballMat->SetPoissonRatio(nu);
	ballMat->SetRestitution(COR);
	ballMat->SetFriction(mu);
#else
	ChSharedPtr<ChMaterialSurface> ballMat;
	ballMat = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	ballMat->SetRestitution(COR);
	ballMat->SetFriction(mu);
#endif

	ChSharedBodyPtr ball(new ChBody(new ChCollisionModelParallel));

	ChVector<> pos(x, y, z);

	ball->SetMaterialSurface(ballMat);

	ball->SetIdentifier(ballId);
	ball->SetMass(ballMass);
	ball->SetInertiaXX(inertia);
	ball->SetPos(pos);
	ball->SetRot(ChQuaternion<>(1, 0, 0, 0));
	ball->SetBodyFixed(false);
	ball->SetCollide(true);

	ball->GetCollisionModel()->ClearModel();
	utils::AddSphereGeometry(ball.get_ptr(), ballRad, ChVector<double>(0, 0, 0), ChQuaternion<double>(1, 0, 0, 0), true);
	ball->GetCollisionModel()->SetFamily(3);
	ball->GetCollisionModel()->BuildModel();

	sys->AddBody(ball);

}

// =============================================================
// III. AddBalls
// Function that calculates position of the sphere in the 3D grid
// After calculating its position, functions adds sphere to the system
// via AddBall(...) function.
// In for (int ii = 1; ii <= zLay; ii++){ ... } hexagonal packing algorithm
// is written. It can be replaced just to produce (x,y,z) for AddBall(...) function

void AddBalls(ChSystemParallel* sys) {

	int ballId = 0;

	double x, y, z;

	for (int ii = 1; ii <= zLay; ii++){
		z = ballRad + (ii - 1) * 2 * sqrt(2) / 3 * ballDiam;

		if (ii % 4 == 0) {
			for (int jj = 1; jj <= yLay - 1; jj++){
				y = -(yLay - 3.0) / 2.0 * sqrt(3) * ballRad - 2. / 3.* sqrt(3) * ballRad + (jj - 1) * sqrt(3) * ballRad;
				if (jj % 2 != 0){
					for (int kk = 1; kk <= xLay - 1; kk++){
						x = (-(xLay - 2) + 2 * (kk - 1)) * ballRad;
						ballId++;
						AddBall(sys, x, y, z, ballId);
					}
				}
				else{
					for (int kk = 1; kk <= xLay - 2; kk++){
						x = (-(xLay - 3) + 2 * (kk - 1)) * ballRad;
						ballId++;
						AddBall(sys, x, y, z, ballId);
					}
				}
			}
		}
		else if (ii % 2 == 0){
			//=====================================
			for (int jj = 1; jj <= yLay - 1; jj++){
				y = -(-(yLay - 3.0) / 2.0 * sqrt(3) * ballRad - 2. / 3.* sqrt(3) * ballRad + (jj - 1) * sqrt(3) * ballRad);
				if (jj % 2 != 0){
					for (int kk = 1; kk <= xLay - 1; kk++){
						x = (-(xLay - 2) + 2 * (kk - 1)) * ballRad;
						ballId++;
						AddBall(sys, x, y, z, ballId);
					}
				}
				else{
					for (int kk = 1; kk <= xLay - 2; kk++){
						x = (-(xLay - 3) + 2 * (kk - 1)) * ballRad;
						ballId++;
						AddBall(sys, x, y, z, ballId);
					}
				}
			}
			//=====================================
		}
		else{
			for (int jj = 1; jj <= yLay; jj++){
				y = -(yLay - 1.0) / 2.0 * sqrt(3) * ballRad + (jj - 1) * sqrt(3) * ballRad;
				if (jj % 2 != 0){
					for (int kk = 1; kk <= xLay; kk++){
						x = (-(xLay - 1) + 2 * (kk - 1)) * ballRad;
						ballId++;
						AddBall(sys, x, y, z, ballId);
					}
				}
				else{
					for (int kk = 1; kk <= xLay - 1; kk++){
						x = (-(xLay - 2) + 2 * (kk - 1)) * ballRad;
						ballId++;
						AddBall(sys, x, y, z, ballId);
					}
				}
			}
		}
	}
}


int main(int argc, char* argv[]) {

	//-----------------------------------------------------------------
	// a) Make output directories
	// b) Create system
	// c) Set the gravity
	// d) Create the OpenGL visualization window
	// e) Set the solver's parameters
	// f) Set the way of collison detection
	// g) Create container
	// h)  Create spheres
	// i) Loop

	// ----------------------------------------------------------------
	// a) Make output directories, output files and show names of the
	// columns the will contain output data
	if (ChFileutils::MakeDirectory(out_dir.c_str()) < 0) {
		cout << "Error creating directory " << out_dir << endl;
		return 1;
	}
	if (ChFileutils::MakeDirectory(pov_dir.c_str()) < 0) {
		cout << "Error creating directory " << pov_dir << endl;
		return 1;
	}

	ChStreamOutAsciiFile forceSensorStream(forceSensor_file.c_str());
	forceSensorStream.SetNumFormat("%16.4e");
	// ----------------------------------------------------------------
	// ----------------------------------------------------------------
	real3 forceSensor;
	// ----------------------------------------------------------------
	// ----------------------------------------------------------------
	cout << "time" << " \t" << "Fx" << "\t" << "Fy" << "\t" << "Fz" << "\n";
	forceSensorStream << "time" << " \t" << "Fx" << "\t" << "Fy" << "\t" << "Fz" << "\n";

	// ----------------------------------------------------------------
	// b) Create the system
#ifdef USE_DEM
	cout << "Create DEM-P system" << endl;
	const std::string title = "Wave propagation in granular media. DEM-P.";
	ChBody::ContactMethod contact_method = ChBody::DEM;
	ChSystemParallelDEM* my_system = new ChSystemParallelDEM();
#else
	cout << "Create DEM-C system" << endl;
	const std::string title = "Wave propagation in granular media. DEM-C.";
	ChBody::ContactMethod contact_method = ChBody::DVI;
	ChSystemParallelDVI* my_system = new ChSystemParallelDVI();
#endif

	// ----------------------------------------------------------------
	// c) Set the gravity
	my_system->Set_G_acc(ChVector<>(0, 0, -gravity));

	// ----------------------------------------------------------------
	// d) Create the OpenGL visualization window
#ifdef CHRONO_PARALLEL_HAS_OPENGL
	opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
	gl_window.Initialize(800, 600, title.c_str(), my_system);
	gl_window.SetCamera(ChVector<>(4 * xDim, 0.001, 3 * zDim), ChVector<>(0, 0, zDim), ChVector<>(0, 0, 1));
	//gl_window.SetCamera(ChVector<>(2 * xDim, 2 * yDim, 4 * zDim), ChVector<>(0, 0, 0), ChVector<>(0, 0, 1));
	gl_window.SetRenderMode(opengl::WIREFRAME);
#endif

	// ----------------------------------------------------------------
	// e) Set the solver's parameters
	int max_threads = my_system->GetParallelThreadNumber();
	if (threads > max_threads) threads = max_threads;
	my_system->SetParallelThreadNumber(threads);
	omp_set_num_threads(threads);

	my_system->GetSettings()->max_threads = threads;
	my_system->GetSettings()->perform_thread_tuning = thread_tuning;

	my_system->GetSettings()->solver.use_full_inertia_tensor = false;
	my_system->GetSettings()->solver.tolerance = tolerance;
	my_system->GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
	my_system->GetSettings()->solver.clamp_bilaterals = clamp_bilaterals;
	my_system->GetSettings()->solver.bilateral_clamp_speed = bilateral_clamp_speed;

#ifdef USE_DEM
	my_system->GetSettings()->solver.contact_force_model = HERTZ;
	my_system->GetSettings()->solver.tangential_displ_mode = ONE_STEP;
#else
	my_system->GetSettings()->solver.solver_mode = SLIDING;
	my_system->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
	my_system->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
	my_system->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
	my_system->GetSettings()->solver.alpha = 0;
	my_system->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
	my_system->ChangeSolverType(APGD);

	// 5% of the radius? No, I think that we dont need it.
	my_system->GetSettings()->collision.collision_envelope = 0.0; 
#endif

	// ----------------------------------------------------------------
	// f) Set the way of collison detection
	my_system->GetSettings()->collision.bins_per_axis = I3(10, 10, 10);
	my_system->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_MPR;

	// ----------------------------------------------------------------
	// g) Create container and add it to the system
	ChSharedPtr<ChBody> wall_1(new ChBody(new ChCollisionModelParallel, contact_method)); // bottom
	//ChSharedPtr<ChBody> wall_2(new ChBody(new ChCollisionModelParallel, contact_method)); // top - optional
	ChSharedPtr<ChBody> wall_3(new ChBody(new ChCollisionModelParallel, contact_method)); // back x side
	ChSharedPtr<ChBody> wall_4(new ChBody(new ChCollisionModelParallel, contact_method)); // front x wall
	ChSharedPtr<ChBody> wall_5(new ChBody(new ChCollisionModelParallel, contact_method)); // back y wall - the one which will be excited
	ChSharedPtr<ChBody> wall_6(new ChBody(new ChCollisionModelParallel, contact_method)); // front y wall
	
	AddWall(my_system, wall_1, -1, ChVector<>(0, 0, -thickness / 2), ChVector<>(xDim, yDim, thickness / 2), true);
	//AddWall(my_system, wall_2, -2, ChVector<>(0, 0, zDim + thickness / 2), ChVector<>(xDim, yDim, thickness / 2), false);
	AddWall(my_system, wall_3, -3, ChVector<>(-(xDim / 2 + thickness / 2), 0, zDim / 2), ChVector<>(thickness / 2, yDim, zDim), false);
	AddWall(my_system, wall_4, -4, ChVector<>(xDim / 2 + thickness / 2, 0, zDim / 2), ChVector<>(thickness / 2, yDim, zDim), false);
	AddWall(my_system, wall_5, -5, ChVector<>(0, -(yDim / 2 + thickness / 2), zDim / 2), ChVector<>(xDim, thickness / 2, zDim), true);
	AddWall(my_system, wall_6, -6, ChVector<>(0, yDim / 2 + thickness / 2, zDim / 2), ChVector<>(xDim, thickness / 2, zDim), false);

	// ----------------------------------------------------------------
	// h) Create spheres and add them to the system
	AddBalls(my_system);

	// ----------------------------------------------------------------
	// i) Loop
	while (my_system->GetChTime() <= simulation_time) {
		

#ifdef VELOCITY_EXCITATION
		//// ________________________________________________________________
		//// Make the wall unfixed and add Link Joint to it
		//if (settled == false && my_system->GetChTime() >= timeStartPushing){
		//	settled = true;
		//	wall_5->SetBodyFixed(false);

		//	// Define quaternion representing:
		//	// - a rotation of -90 degrees around x (z2y)
		//	ChQuaternion<> z2y;
		//	z2y.Q_from_AngAxis(-CH_C_PI / 2, ChVector<>(1, 0, 0));

		//	ChSharedPtr<ChLinkLockPrismatic> prismatic_wall_5_1(new ChLinkLockPrismatic);
		//	prismatic_wall_5_1->SetName("prismatic_wall_5_1");
		//	prismatic_wall_5_1->Initialize(wall_5, wall_1, ChCoordsys<>(ChVector<>(0, 0, 0), z2y));
		//	my_system->AddLink(prismatic_wall_5_1);
		//}

		// ________________________________________________________________
		// Activate wall to start pushing the specimen
		// After timePushing make the wall comming back
		// After whole the action wall becomes fixed again.
		// ITS ONLY ONE OF THE WAYS TO EXCITE THE SPECIMEN WITH VELOCITY. REMAINING WAYS:
		// 1. Uncomment: 
		//		a) "Make the wall unfixed and add Link Joint to it"
		//		b) wall_5->SetPos_dt(ChVector<>(0, velocity, 0));
		//		c) wall_5->SetBodyFixed(true);
		// But that way has a drawback: wall is unfixed and spheres deccelerate it, thus applied velocity is smaller than expected
		// 2. ...

		if (my_system->GetChTime() <= timeStartPushing + timePushing && my_system->GetChTime() >= timeStartPushing){
			wall_5->SetPos(wall_5->GetPos() + ChVector<>(0, velocity * time_step, 0));
			//wall_5->SetPos_dt(ChVector<>(0, velocity, 0));
		}
		else if ((my_system->GetChTime() > timeStartPushing + timePushing) 
			&& (my_system->GetChTime() <= timeStartPushing + 2 * timePushing)){
			wall_5->SetPos(wall_5->GetPos() - ChVector<>(0, velocity * time_step, 0));
			//wall_5->SetPos_dt(ChVector<>(0, -velocity, 0));
		}
		else{
			//wall_5->SetBodyFixed(true);
		}
#else
		// ________________________________________________________________
		// Make the wall unfixed and add Link Joint to it
		if (settled == false && my_system->GetChTime() >= timeStartPushing){
			settled = true;
			wall_5->SetBodyFixed(false);

			// Define quaternion representing:
			// - a rotation of -90 degrees around x (z2y)
			ChQuaternion<> z2y;
			z2y.Q_from_AngAxis(-CH_C_PI / 2, ChVector<>(1, 0, 0));

			ChSharedPtr<ChLinkLockPrismatic> prismatic_wall_5_1(new ChLinkLockPrismatic);
			prismatic_wall_5_1->SetName("prismatic_wall_5_1");
			prismatic_wall_5_1->Initialize(wall_5, wall_1, ChCoordsys<>(ChVector<>(0, 0, 0), z2y));
			my_system->AddLink(prismatic_wall_5_1);
		}

		if (my_system->GetChTime() <= timeStartPushing + timePushing && my_system->GetChTime() >= timeStartPushing){
			wall_5->Empty_forces_accumulators();
			wall_5->Accumulate_force(ChVector<>(0, force, 0), wall_5->GetPos(), false);
		}
		if (pushed == false && my_system->GetChTime() >= timeStartPushing + 2*timePushing){
			pushed = true;
			wall_5->SetBodyFixed(true);
		}


#endif

		// ________________________________________________________________
		//  Do time step
#ifdef CHRONO_PARALLEL_HAS_OPENGL
		if (gl_window.Active()) {
			gl_window.DoStepDynamics(time_step);
			gl_window.Render();
		}
		else
			break;
#else
		my_system->DoStepDynamics(time_step);
#endif

		// ________________________________________________________________
		// Produce output data
		if (my_system->GetChTime() >= data_out_frame * data_out_step) {

#ifndef USE_DEM
			my_system->CalculateContactForces();
#endif
			uint a = wall_6->GetId();
			forceSensor = my_system->GetBodyContactForce(a);
			cout << my_system->GetChTime() << "\t" << forceSensor.x << "\t" << forceSensor.y << "\t" << forceSensor.z << "\n";
			forceSensorStream << my_system->GetChTime() << "\t" << forceSensor.x << "\t" << forceSensor.y << "\t" << forceSensor.z << "\n";
			data_out_frame++;
		}

		// ________________________________________________________________
		// Produce POV-Ray data
		if (write_povray_data && my_system->GetChTime() >= visual_out_frame * visual_out_step) {
			char filename[100];
			sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), visual_out_frame + 1);
			utils::WriteShapesPovray(my_system, filename, false);

			visual_out_frame++;
		}

	}// end while

	return 0;
}
