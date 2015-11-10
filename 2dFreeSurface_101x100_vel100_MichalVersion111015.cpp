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

#include "unit_PARALLEL/chrono_utils/ChUtilsCreators.h"
#include "unit_PARALLEL/chrono_utils/ChUtilsGenerators.h"
#include "unit_PARALLEL/chrono_utils/ChUtilsInputOutput.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::flush;
using std::endl;

// =============================================================
// Above headers and namespaces
#define HEXAGONAL_CLOSE_2D

#define VELOCITY_EXCITATION

#define USE_DEM

#undef CHRONO_OPENGL

#ifdef CHRONO_OPENGL
	#include "chrono_opengl/ChOpenGLWindow.h"
#endif

double geomScale = 1e6; //1;

// Parameters you set up to have different simulations
double overlap = 1e-3;//1e-9 * geomScale; //
double time_step = 1e-13;//1.0e-13;
double velocity = 100 * geomScale;// excitation velocity || m/s

double simulation_time = 5e-8;
bool write_povray_data = true;
double data_out_step = 50 * time_step;//1e-8;       // time interval between data outputs
double visual_out_step = 50 * time_step;//1e-7;     // time interval between PovRay outputs
double energy_out_step = 1000 * time_step; // interval between output conderning potential and kinetoc energy

#ifdef USE_DEM
	const std::string out_dir1 = "Seattle_3";
	const std::string out_dir = out_dir1 + "/2dFreeSurface_101x100_vel=100_noGrav";
#endif

const std::string pov_dir = out_dir + "/POVRAY";

const std::string pov_vel_dir = out_dir + "/POVRAY_velocities";
const std::string energy_dir = out_dir + "/Energy";
const std::string particles_dir = out_dir + "/Particles";

const std::string time_file = out_dir + "/time.dat";

const std::string checkPoint_file = out_dir + "/checkPoint.dat";

void OutputData(ChSystemParallel* sys, int out_frame, double time) {
	char filename[100];
	sprintf(filename, "%s/data_%04d.dat", pov_dir.c_str(), out_frame);
	utils::WriteShapesPovray(sys, filename, false);
}

double gravity = 0;//9.81 * geomScale; // m/s/s

double ballRad = 0.5;//0.5e-6 * geomScale; // micrometer
double ballDiam = 2.0 * ballRad;

float mu = 0.18f; // needless in Seattle model
float ballDensity = 2000.0f / geomScale / geomScale / geomScale;//1.0f; // kg/m/m/m
float ballMass = ballDensity * 4.0f / 3.0f * CH_C_PI * ballRad * ballRad * ballRad;

float E = 73e9 / geomScale; //1e9 / geomScale;//1e9;// / Young's modulus. Pa = N/m/m = kg/s/s/m || to watch out for crazy vibrations
float nu = 0.17f;//0.32f;

ChVector<> inertia = (2.0 / 5.0) * ballMass * ballRad * ballRad * ChVector<>(1, 1, 1);

#ifdef HEXAGONAL_CLOSE_2D
	double xLay = 101;
	double zLay = 100;
	double yLay = 1;

	double pushedBall = (xLay-1) / 2;
#endif

// xOverlap - for 2D and 3D only
// yOverlap above - different values for ONE_D_ARRANGEMENT and different for HEXAGONAL_CLOSE because:
	// in ONE_D_ARRANGEMENT - spheres touch each other in the upright way
	// in HEXAGONAL_CLOSE - they touch each other in the "triangle" way (equilateral triangle to be clear)
// zOverlap - for 3D only

//	ChVector<> CameraLocation;
//	ChVector<> CameraLookAt;

#ifdef HEXAGONAL_CLOSE_2D
	// add some overlap to here
	double xDim = (ballDiam - overlap) * xLay + 0.5 *(ballDiam - overlap);
	double yDim = ballDiam;
	double zDim = (zLay - 1) * (ballDiam - overlap) * sqrt(3) / 2;

	ChVector<> CameraLocation = ChVector<>(xDim / 2.0, 10 * yDim, 0);
	ChVector<> CameraLookAt = ChVector<>(xDim / 2.0, 0.01, 0.01);
#endif

int threads = 2;

bool thread_tuning = false;

#ifdef USE_DEM
	//double time_step = 1.0e-9; // is was set up above
	double tolerance = 1.0;//0.001;   // Arman, check this
	int max_iteration_bilateral = 50;
#endif

bool clamp_bilaterals = false;
double bilateral_clamp_speed = 0.1;

double timeStartPushing = 0;//50*time_step;//3.0e-6; // start time for pushing the wall

int data_out_frame = 0;
int visual_out_frame = 0;
int energy_out_frame = 0;

//#ifndef ONE_D_ARRANGEMENT
//	void AddWall(ChSystemParallel* sys, ChSharedPtr<ChBody> wall,
//		int wallId, ChVector<> pos, ChVector<> dim,
//		bool ifFixed, bool visualization){
//
//		wall->SetIdentifier(wallId);
//		//double wallMass = dim.x * dim.y * dim.z * wallDensity;
//		wall->SetMass(mass_ext);
//		wall->SetPos(pos);
//		wall->SetBodyFixed(ifFixed);
//		wall->SetCollide(true);
//
//		#ifdef USE_DEM
//			ChSharedPtr<ChMaterialSurfaceDEM> mat_ext;
//			mat_ext = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
//			mat_ext->SetYoungModulus(Y_ext);
//			mat_ext->SetPoissonRatio(nu_ext);
//			mat_ext->SetRestitution(COR_ext);
//			mat_ext->SetFriction(mu_ext);
//		#else
//			ChSharedPtr<ChMaterialSurface> mat_ext;
//			mat_ext = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
//			mat_ext->SetRestitution(COR_ext);
//			mat_ext->SetFriction(mu_ext);
//		#endif
//
//		wall->SetMaterialSurface(mat_ext);
//
//		wall->GetCollisionModel()->ClearModel();
//
//		utils::AddBoxGeometry(wall.get_ptr(), dim, ChVector<>(0, 0, 0), ChQuaternion<double>(1, 0, 0, 0), visualization);
//		wall->GetCollisionModel()->SetFamily(1);
//		wall->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
//
//		wall->GetCollisionModel()->BuildModel();
//		sys->AddBody(wall);
//	}
//#endif

void AddBall(ChSystemParallel* sys, double x, double y, double z, int ballId, bool ifFixed){

	#ifdef USE_DEM
		ChSharedPtr<ChMaterialSurfaceDEM> ballMat;
		ballMat = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
		ballMat->SetYoungModulus(E);
		ballMat->SetPoissonRatio(nu);
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
	ball->SetBodyFixed(ifFixed);
	ball->SetCollide(true);

	ball->GetCollisionModel()->ClearModel();
	utils::AddSphereGeometry(ball.get_ptr(), ballRad, ChVector<double>(0, 0, 0), ChQuaternion<double>(1, 0, 0, 0), true);
	ball->GetCollisionModel()->SetFamily(1);
	ball->GetCollisionModel()->SetDefaultSuggestedEnvelope(0.00);
	ball->GetCollisionModel()->SetDefaultSuggestedMargin(0.00 * ballDiam);


	ball->GetCollisionModel()->BuildModel();

	sys->AddBody(ball);
}

#ifdef HEXAGONAL_CLOSE_2D

void AddBalls(ChSystemParallel* sys) {

	int ballId = 0;

	bool ifFixed = false; // before figuring out periodic boundaries let spheres at the boundaries be fixed

	double x = 0., z = 0.;
	double y = 0.;
	double dimention = ballDiam - overlap; // dimention = length of the triangle's side (it repeats all the time in the code below - that is why it is no sense in computing it all the time)

	for (int ii = 0; ii < zLay; ii++){

		z = ii * (- dimention) * sqrt(3) / 2; // height of the triangle

		for (int jj = 0; jj < xLay; jj++){

			// ========= if Fixed ??? =================
			ifFixed = false;

			if (ii == zLay - 1) // last row should be fixed
				ifFixed = true;

			if ((jj == 0) || (jj == xLay - 1)) // first and last column should be fixed
				ifFixed = true;
			// =========================================

			if (ii % 2 == 0){
				x = jj * dimention;
				ballId++;
				AddBall(sys, x, y, z, ballId, ifFixed);
			}
			else{
				x = jj * dimention + dimention * 0.5;
				ballId++;
				AddBall(sys, x, y, z, ballId, ifFixed);
			}
		}
	}
}

#endif

int main(int argc, char* argv[]){

	// ================================================
	// ================================================
	
	const int howManyParticles = 6;
	int particleNumber[howManyParticles];
	particleNumber[0] = 0; // the middle one (0th) in the first row (that is why there is an odd amount of spheres in the first row)
	particleNumber[1] = 10;
	particleNumber[2] = 20;
	particleNumber[3] = 30;
	particleNumber[4] = 40;
	particleNumber[5] = 50;

	// recalculating particles' numbers for their real ID
	const int numOfParticles = 2 * howManyParticles - 1; // -1 because particle 0 is the same particle in row 0... we don't want to have same output twice
	int particleID[numOfParticles];
	for (int ii = 0; ii < howManyParticles; ii++){
		particleID[ii] = (xLay - 1) / 2 + particleNumber[ii];
	}
	for (int ii = 0; ii < howManyParticles - 1; ii++){ // note that number of iterations in the loop changed
		particleID[ii + howManyParticles] = (particleNumber[ii + 1] - 1)*xLay + (xLay - 1) / 2;
	}

	// ================================================
	// ================================================

	//cout << "ballDensity" << ballDensity;

	if (ChFileutils::MakeDirectory(out_dir1.c_str()) < 0) {
		cout << "Error creating directory " << out_dir1 << endl;
		return 1;
	}
	if (ChFileutils::MakeDirectory(out_dir.c_str()) < 0) {
		cout << "Error creating directory " << out_dir << endl;
		return 1;
	}
	if (ChFileutils::MakeDirectory(pov_dir.c_str()) < 0) {
		cout << "Error creating directory " << pov_dir << endl;
		return 1;
	}
	if (ChFileutils::MakeDirectory(pov_vel_dir.c_str()) < 0) {
		cout << "Error creating directory " << pov_vel_dir << endl;
		return 1;
	}
	if (ChFileutils::MakeDirectory(energy_dir.c_str()) < 0) {
		cout << "Error creating directory " << energy_dir << endl;
		return 1;
	}
	if (ChFileutils::MakeDirectory(particles_dir.c_str()) < 0) {
		cout << "Error creating directory " << particles_dir << endl;
		return 1;
	}

#ifdef USE_DEM
	cout << "Create DEM-P system" << endl;
	const std::string title = "Wave propagation in granular media. DEM-P.";
	ChSystemParallelDEM* my_system = new ChSystemParallelDEM();
#endif

	my_system->Set_G_acc(ChVector<>(0, 0, -gravity));

#ifdef CHRONO_OPENGL
	opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
	gl_window.Initialize(800, 600, title.c_str(), my_system);

	gl_window.SetCamera(CameraLocation, ChVector<>(0.2, 0.3, 0.4), ChVector<>(0, 0, 1), ballDiam, ballDiam);
//	gl_window.viewer->render_camera.camera_scale = 2.0/(1000.0);
//	gl_window.viewer->render_camera.near_clip = .001;
	//gl_window.SetRenderMode(opengl::WIREFRAME);

#endif

	int max_threads = my_system->GetParallelThreadNumber();
	if (threads > max_threads)
		threads = max_threads;
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
	my_system->GetSettings()->solver.contact_force_model = SEATTLE;
	my_system->GetSettings()->solver.tangential_displ_mode = MULTI_STEP;
#endif

	my_system->GetSettings()->collision.bins_per_axis = I3(20, 1, 20);
	my_system->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_MPR;
	
	AddBalls(my_system);

	bool pushed = false;

	while (my_system->GetChTime() <= simulation_time) {

#ifdef VELOCITY_EXCITATION
		//// push the first sphere after the certain time
		if (pushed == false && my_system->GetChTime() >= timeStartPushing){
		
			pushed = true;
			ChSharedPtr<ChBody> ball_X;
			ball_X = ChSharedPtr<ChBody>(my_system->Get_bodylist()->at(0 + pushedBall));
			my_system->Get_bodylist()->at(0 + pushedBall)->AddRef();
			ball_X->SetPos_dt(ball_X->GetPos_dt() + ChVector<>(0, 0, - velocity ));
			//ball_X->SetPos_dt(ball_X->GetPos() + ChVector<>(0, velocity*time_step, 0));
		}
#else
		//Force excitation
#endif

#ifdef CHRONO_OPENGL
		if (gl_window.Active()) {
			gl_window.DoStepDynamics(time_step);
			gl_window.Render();
		}
		else{
			break;
		}
#else
		my_system->DoStepDynamics(time_step);
#endif


#ifdef HEXAGONAL_CLOSE_2D
		{ // projections + information about violation
			ChSharedPtr<ChBody> ball_X;

			// .. and some checking of plane violation
			double planeViolation = 0;
			double planeViolation2 = 0;
			// 
			double planeViolationAfter = 0;
			double planeViolationAfter2 = 0;


			for (int ii = 0; ii < xLay*zLay; ii++){
				ball_X = ChSharedPtr<ChBody>(my_system->Get_bodylist()->at(ii));
				my_system->Get_bodylist()->at(ii)->AddRef();

				planeViolation += sqrt(ball_X->GetPos().y*ball_X->GetPos().y);
				planeViolation2 += ball_X->GetPos().y;


				ball_X->SetPos( ChVector<double>( ball_X->GetPos().x, 0.0, ball_X->GetPos().z ) );
				ball_X->SetRot(ChQuaternion<double>(ball_X->GetRot().e0, 0.0, ball_X->GetRot().e2, 0.0));
				ball_X->SetPos_dt(ChVector<double>(ball_X->GetPos_dt().x, 0.0, ball_X->GetPos_dt().z));
				ball_X->SetRot_dt(ChQuaternion<double>( - 0.5 * ball_X->GetRot().e2 * ball_X->GetWvel_loc().y,
					0.0,
					0.5 * ball_X->GetRot().e2 * ball_X->GetWvel_loc().y,
					0.0));

				planeViolationAfter += sqrt(ball_X->GetPos().y*ball_X->GetPos().y);
				planeViolationAfter2 += ball_X->GetPos().y;


			}

			char planeFile[100];
			sprintf(planeFile, "%s/planeViolation.dat", out_dir.c_str());

			FILE * planeViol = fopen(planeFile, "a");
			fprintf(planeViol, "%.14f,%.14f,%.14f,%.14f,\n", planeViolation, planeViolation2, planeViolationAfter, planeViolationAfter2);
			fclose(planeViol);
		}
#endif

		if (my_system->GetChTime() >= energy_out_frame * energy_out_step) {

			char energyFileName[100];
			sprintf(energyFileName, "%s/energy%04d.dat", energy_dir.c_str(), energy_out_frame + 1);
			FILE * energyFile = fopen(energyFileName, "w");

			int licznik = 0;
			for (int jj = 0; jj < xLay; jj++){
				for (int ii = 0; ii < zLay; ii++){

					ChSharedPtr<ChBody> ball_X;
					ball_X = ChSharedPtr<ChBody>(my_system->Get_bodylist()->at(licznik));
					my_system->Get_bodylist()->at(licznik)->AddRef();

					double height = ball_X->GetPos().z;
					double velMag = sqrt(ball_X->GetPos_dt().x*ball_X->GetPos_dt().x + ball_X->GetPos_dt().z *ball_X->GetPos_dt().z);
					double angVel = ball_X->GetWvel_par().y;

					fprintf(energyFile, "%.13f,%.13f,%.13f,", height, velMag, angVel);
					licznik++;
				}
				fprintf(energyFile, "\n");
			}
			fclose(energyFile);
			energy_out_frame++;
		}


		if (my_system->GetChTime() >= data_out_frame * data_out_step) {

#ifndef USE_DEM
			my_system->CalculateContactForces();
#endif

			// _________________________________________________________________________
			// time information
			FILE * timeWrite = fopen(time_file.c_str(), "a");
			fprintf(timeWrite, "%.13f,\n", my_system->GetChTime());
			fclose(timeWrite);

			// _________________________________________________________________________
			// information about certain particles' position and velocities
		
			for (int jj = 0; jj < numOfParticles; jj++){

				char rowFileName[100];
				sprintf(rowFileName, "%s/particle%04d.dat", particles_dir.c_str(), particleID[jj]);
				FILE * rowFile = fopen(rowFileName, "a");

				ChSharedPtr<ChBody> ball_X;
				ball_X = ChSharedPtr<ChBody>(my_system->Get_bodylist()->at(particleID[jj]));
				my_system->Get_bodylist()->at(particleID[jj])->AddRef();

				fprintf(rowFile, "%.13f,%.13f,%.13f,%.13f,%.13f,%.13f,%.13f,\n",
					my_system->GetChTime(),
					ball_X->GetPos().x, ball_X->GetPos().z, ball_X->GetRotAngle(),
					ball_X->GetPos_dt().x, ball_X->GetPos_dt().z, ball_X->GetWvel_par().y);

				fclose(rowFile);
			}
			data_out_frame++;
		} // end of if - writing data to file

		// ________________________________________________________________
		// Produce POV-Ray data
		if (write_povray_data && my_system->GetChTime() >= visual_out_frame * visual_out_step) {


			// =====================================================================
			// =====================================================================

			char filename[100];
			sprintf(filename, "%s/data_%04d.dat", pov_dir.c_str(), visual_out_frame + 1);
			utils::WriteShapesPovray(my_system, filename, false);

			// =====================================================================
			// =====================================================================

			// ... and velocity magnitude of spheres
			char AllVel_file[100];
			sprintf(AllVel_file, "%s/data_%04d.dat", pov_vel_dir.c_str(), visual_out_frame + 1);
			FILE * AllVel_a = fopen(AllVel_file, "w");

			double maxVelocity = 0;

			ChSharedPtr<ChBody> ball_X;
			for (int ii = 0; ii < xLay*zLay; ii++){
				ball_X = ChSharedPtr<ChBody>(my_system->Get_bodylist()->at(ii));
				my_system->Get_bodylist()->at(ii)->AddRef();

				double velo = sqrt(ball_X->GetPos_dt().x*ball_X->GetPos_dt().x +
					ball_X->GetPos_dt().y*ball_X->GetPos_dt().y +
					ball_X->GetPos_dt().z*ball_X->GetPos_dt().z);

				fprintf(AllVel_a, "%f,\n", velo);

				if (maxVelocity < velo){
					maxVelocity = velo;
				}
			}

			fclose(AllVel_a);

			// =====================================================================
			// =====================================================================

			// maximal velocity
			char MaxVelocity_file[100];
			sprintf(MaxVelocity_file, "%s/maxVelocity_%04d.dat", pov_vel_dir.c_str(), visual_out_frame + 1);
			FILE * MaxVelocity_a = fopen(MaxVelocity_file, "w");
			fprintf(MaxVelocity_a, "%f,\n", maxVelocity);
			fclose(MaxVelocity_a);

			// =====================================================================
			// =====================================================================

			// ... colors - just to have sth to color them with
			char AllVelColor_file[100];
			sprintf(AllVelColor_file, "%s/dataColor_%04d.dat", pov_vel_dir.c_str(), visual_out_frame + 1);
			FILE * AllVelColor_a = fopen(AllVelColor_file, "w");

			fprintf(AllVelColor_a, "%.9f,\n",
				sqrt(ball_X->GetPos_dt().x*ball_X->GetPos_dt().x +
				ball_X->GetPos_dt().y*ball_X->GetPos_dt().y +
				ball_X->GetPos_dt().z*ball_X->GetPos_dt().z) / velocity);

			for (int ii = 0; ii < xLay*zLay; ii++){
				ball_X = ChSharedPtr<ChBody>(my_system->Get_bodylist()->at(ii));
				my_system->Get_bodylist()->at(ii)->AddRef();

				double velo = sqrt(ball_X->GetPos_dt().x*ball_X->GetPos_dt().x +
					ball_X->GetPos_dt().y*ball_X->GetPos_dt().y +
					ball_X->GetPos_dt().z*ball_X->GetPos_dt().z);

				fprintf(AllVelColor_a, "%.9f,\n", velo / velocity);
			}

			fclose(AllVelColor_a);

			visual_out_frame++;
		}
	}// end of while

	cout << "Write checkpoint data to " << flush;
	utils::WriteCheckpoint(my_system, checkPoint_file);

	return 0;
	
}