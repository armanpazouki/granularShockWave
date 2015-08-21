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

#define m2mm 1e6//1e6//1e6//1//1000000
#define m2mm2 1//1e-6//1//1 

//#define HEXAGONAL_CLOSE
#define LINK

#define VELOCITY_EXCITATION

#define USE_DEM

//#undef CHRONO_PARALLEL_HAS_OPENGL

#ifdef CHRONO_PARALLEL_HAS_OPENGL
	#include "chrono_opengl/ChOpenGLWindow.h"
#endif

// Parameters you set up to have different simulations
double yOverlap = 0.001 * m2mm2; // 
double time_step = 1e-13;//1.0e-13;
double velocity = 5 * m2mm;// * m2mm; // excitation velocity || m/s

double simulation_time = 2e-9;
bool write_povray_data = false;
double data_out_step = 10*time_step;//1e-8;       // time interval between data outputs
double visual_out_step = 100*time_step;//1e-7;     // time interval between PovRay outputs

#ifdef USE_DEM
	//const std::string out_dir = "../Seattle_adhesion_DEM_ov=0.1_ts=1.0e-9_tp=60ts_vel=25_out=1.0e-8";
	const std::string out_dir = "../FinalOnes2_noScaling_withMinusByinf/vel=5";
	//const std::string out_dir = "../Seattle_adhesion_DEM_test3week";
#else
	const std::string out_dir = "../Seattle_adhesion_DVI";
#endif

const std::string pov_dir = out_dir + "/POVRAY";

//const std::string particlePos1_file = out_dir + "/particlePosVel1.dat";
const std::string particlePos2_file = out_dir + "/particlePosVel2.dat";
const std::string particlePos3_file = out_dir + "/particlePosVel3.dat";

const std::string allPos_file = out_dir + "/allPosVel.dat";
//const std::string particlePos3_file = out_dir + "/particlePosVel3.dat";
const std::string time_file = out_dir + "/time.dat";

const std::string forceSensor_file = out_dir + "/forceSensor.dat";
const std::string checkPoint_file = out_dir + "/checkPoint.dat";

void OutputData(ChSystemParallel* sys, int out_frame, double time) {
	char filename[100];
	sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame);
	utils::WriteShapesPovray(sys, filename, false);
	//cout << "time = " << time << flush << endl;
}

double gravity = 0.0 * m2mm;//9.81 * m2mm; // m/s/s

double ballRad = 0.5 * m2mm2; // micrometer
double ballDiam = 2 * ballRad;

float mu = 0.18f; // needless in Seattle model
float COR = 0.87f; // needless in Seattle model
float ballDensity = 1.0f / m2mm / m2mm / m2mm; // kg/m/m/m
float ballMass = ballDensity * 4.0f / 3.0f * CH_C_PI * ballRad * ballRad * ballRad;

float Y = 1e9;// / m2mm; // Pa = N/m/m = kg/s/s/m || to watch out for crazy vibrations - Y is 1000 times smaller
float nu = 0.32f;

ChVector<> inertia = (2.0 / 5.0) * ballMass * ballRad * ballRad * ChVector<>(1, 1, 1);

#ifdef LINK
	double xLay = 1;
	double zLay = 1;
	double yLay = 22; // the lenght of the link, others have to be = 1.0

	int firstBall = 1 - 1;
	//int middleBall = (int)(yLay + 1) / 2 - 1;
	int lastBall = (int)yLay - 1;
	int fifthBall = 5; // do not substract 1 because the first ball (index 0) is the wall!!!
	int secondBall = 1;
#endif
#ifdef HEXAGONAL_CLOSE
	double xLay = 10;
	double zLay = 25;
	double yLay = 1;
#endif

// xOverlap - for 2D and 3D only
// yOverlap above - different values for LINK and different for HEXAGONAL_CLOSE because:
	// in LINK - spheres touch each other in the upright way
	// in HEXAGONAL_CLOSE - they touch each other in the "triangle" way (equilateral triangle to be clear)
// zOverlap - for 3D only

#ifdef LINK
	double xDim = ballDiam * xLay; // xLay = 1.0 so no need to multiply it by that xLay
	double yDim = yLay * ballDiam - (yLay - 1) * yOverlap; // overlap
	double zDim = ballDiam + (zLay - 1.0) * 2.0 * sqrt(2.0) / 3.0 * ballDiam; // zLay = 1.0 so second component is unnecessary
#endif

#ifdef HEXAGONAL_CLOSE
	// add some overlap to here
	double xDim = ballDiam * xLay;
	double yDim = yLay * ballDiam - (yLay - 1) * ballDiam * yOverlap;
	double zDim = ballDiam + (zLay - 1.0) * 2.0 * sqrt(2.0) / 3.0 * ballDiam;
#endif

//float mu_ext = 0.18f;
//float COR_ext = 0.87f;
//float density_ext = 1060.0f / m2mm / m2mm / m2mm; // kg/m/m/m
//float mass_ext = ballMass;
//float Y_ext = 4.04e9 / m2mm; // Pa = N/m/m = kg/s/s/m || to watch out for crazy vibrations - Y is 1000 times smaller
//float nu_ext = 0.32f;

int threads = 4;

bool thread_tuning = false;

#ifdef USE_DEM
	//double time_step = 1.0e-9;
	double tolerance = 0.1;
	int max_iteration_bilateral = 1000;
#else
	double time_step = 1e-6;
	double tolerance = 0.1;
	int max_iteration_normal = 0;
	int max_iteration_sliding = 1000;
	int max_iteration_spinning = 0;
	int max_iteration_bilateral = 100;
	double contact_recovery_speed = 10e30;
#endif

bool clamp_bilaterals = false;
double bilateral_clamp_speed = 0.1;

double timeStartPushing = 0;//50*time_step;//3.0e-6; // start time for pushing the wall
#ifdef VELOCITY_EXCITATION
	#ifdef USE_DEM
		//double velocity = 100 * m2mm; // velocity of the first wall || m/s // it's declared above (at the beginning of the code)
	#else
		double velocity = 1.5 / 10 * m2mm; // velocity of wall
	#endif
#else
	// Force excitation to be think up
#endif

//bool write_povray_data = true;
//
//double data_out_step = 1e-7;       // time interval between data outputs
//double visual_out_step = 1e-7;     // time interval between PovRay outputs

int data_out_frame = 0;
int visual_out_frame = 0;

#ifndef LINK
	void AddWall(ChSystemParallel* sys, ChSharedPtr<ChBody> wall,
		int wallId, ChVector<> pos, ChVector<> dim,
		bool ifFixed, bool visualization){

		wall->SetIdentifier(wallId);
		//double wallMass = dim.x * dim.y * dim.z * wallDensity;
		wall->SetMass(mass_ext);
		wall->SetPos(pos);
		wall->SetBodyFixed(ifFixed);
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

		utils::AddBoxGeometry(wall.get_ptr(), dim, ChVector<>(0, 0, 0), ChQuaternion<double>(1, 0, 0, 0), visualization);
		wall->GetCollisionModel()->SetFamily(1);
		wall->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);

		wall->GetCollisionModel()->BuildModel();
		sys->AddBody(wall);
	}
#endif

void AddBall(ChSystemParallel* sys, double x, double y, double z, int ballId){

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

#ifdef HEXAGONAL_CLOSE
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
#endif

#ifdef LINK
	void AddBalls(ChSystemParallel* sys){

		int ballId = 0;

		double x, y, z;
		x = 0.0;
		z = ballRad;
		for (int jj = 1; jj <= yLay; jj++){
			//y = -(yLay - 1) * ballRad + (yLay - 1) / 2 * yOverlap +  (jj - 1) * ballDiam - (jj - 1) * yOverlap; // overlap
			y = (jj - 1)*ballDiam - (jj - 1)*yOverlap;
			ballId++;
			AddBall(sys, x, y, z, ballId);
		}
	}
#endif

int main(int argc, char* argv[]){

	cout << "ballDensity" << ballDensity;

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

	real3 forceSensor5;
	real3 forceSensor6;

	//cout << "time" << " \t" << "Fx" << "\t" << "Fy" << "\t" << "Fz" << "\n";
	forceSensorStream << "time" << " \t" << "Fx" << "\t" << "Fy" << "\t" << "Fz" << "\n";

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

	my_system->Set_G_acc(ChVector<>(0, 0, -gravity));

#ifdef CHRONO_PARALLEL_HAS_OPENGL
	opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
	gl_window.Initialize(800, 600, title.c_str(), my_system);

	gl_window.SetCamera(ChVector<>(0, 0, 18 * zDim), ChVector<>(xDim / 1000, 0, 0), ChVector<>(0, 0, 1)); // top view
	//gl_window.SetCamera(ChVector<>(20 * zDim, 0, 0), ChVector<>(0, 0, xDim / 1000 ), ChVector<>(0, 0, 1));// side view
	//gl_window.SetCamera(ChVector<>(- 4 * xDim, 0.001, 3 * zDim), ChVector<>(0, 0, zDim), ChVector<>(0, 0, 1));
	//gl_window.SetCamera(ChVector<>(2 * xDim, 2 * yDim, 4 * zDim), ChVector<>(0, 0, 0), ChVector<>(0, 0, 1));
	//gl_window.SetCamera(ChVector<>(8* xDim, 0.001, 3 * zDim), ChVector<>(0, 0, zDim), ChVector<>(0, 0, 1));
	//gl_window.SetCamera(ChVector<>(2 * xDim, 2 * yDim, 4 * zDim), ChVector<>(0, 0, 0), ChVector<>(0, 0, 1));
	gl_window.SetRenderMode(opengl::WIREFRAME);
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

	my_system->GetSettings()->collision.bins_per_axis = I3(1, 4, 1);
	my_system->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_MPR;

	// the first and the last sphere is fixed
	AddBalls(my_system);
	ChSharedPtr<ChBody> ball_1;
	ball_1 = ChSharedPtr<ChBody>(my_system->Get_bodylist()->at(0 + firstBall)); // 0  because no wall is in there
	my_system->Get_bodylist()->at(0 + firstBall)->AddRef();
	ball_1->SetBodyFixed(true);

	ChSharedPtr<ChBody> ball_N;
	ball_N = ChSharedPtr<ChBody>(my_system->Get_bodylist()->at(0 + lastBall)); // 0  because no wall is in there
	my_system->Get_bodylist()->at(0 + lastBall)->AddRef();
	ball_N->SetBodyFixed(true);

	bool pushed = false;

	//int jj = 1;
	while (my_system->GetChTime() <= simulation_time) {

#ifdef VELOCITY_EXCITATION
		// push the first sphere after the certain time
		if (pushed == false && my_system->GetChTime() >= timeStartPushing){
		
			pushed = true;
			ChSharedPtr<ChBody> ball_X;
			ball_X = ChSharedPtr<ChBody>(my_system->Get_bodylist()->at(0 + firstBall + 1));
			my_system->Get_bodylist()->at(0 + firstBall + 1)->AddRef();
			ball_X->SetPos_dt(ball_X->GetPos_dt() + ChVector<>(0, velocity, 0));
			//ball_X->SetPos_dt(ball_X->GetPos() + ChVector<>(0, velocity*time_step, 0));
		}
#else
		//Force excitation
#endif

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


		if (my_system->GetChTime() >= data_out_frame * data_out_step) {

#ifndef USE_DEM
			my_system->CalculateContactForces();
#endif
			uint a = ball_N->GetId();
			forceSensor6 = my_system->GetBodyContactForce(a);

			uint b = ball_1->GetId();
			forceSensor5 = my_system->GetBodyContactForce(b);

			//cout << my_system->GetChTime() << "\t" << forceSensor6.x << "\t" << forceSensor6.y << "\t" << forceSensor6.z << "\t"
			//	<< forceSensor5.x << "\t" << forceSensor5.y << "\t" << forceSensor5.z << "\n";
			forceSensorStream << my_system->GetChTime() << "\t" << forceSensor6.x << "\t" << forceSensor6.y << "\t" << forceSensor6.z << "\t"
				<< forceSensor5.x << "\t" << forceSensor5.y << "\t" << forceSensor5.z << "\n";


			// _________________________________________________________________________
			// time information

			FILE * timeWrite = fopen(time_file.c_str(), "a");
			fprintf(timeWrite, "%.14f,\n", my_system->GetChTime());
			fclose(timeWrite);

			//// _________________________________________________________________________
			//// first row information

			//FILE * PosVelWrite1 = fopen(particlePos1_file.c_str(), "a");
			//{
			//	ChSharedPtr<ChBody> ball_X;
			//	ball_X = ChSharedPtr<ChBody>(my_system->Get_bodylist()->at(0 + firstBall + 1));
			//	my_system->Get_bodylist()->at(0 + firstBall + 1)->AddRef();

			//	fprintf(PosVelWrite1, "%f,%f,%f,%f,%f,%f,\n",
			//		ball_X->GetPos().x, ball_X->GetPos().y, ball_X->GetPos().z,
			//		ball_X->GetPos_dt().x, ball_X->GetPos_dt().y, ball_X->GetPos_dt().z);
			//}

			//fclose(PosVelWrite1);
			//// _________________________________________________________________________
			//// middle row information

			//FILE * PosVelWrite2 = fopen(particlePos2_file.c_str(), "a");
			//{
			//	ChSharedPtr<ChBody> ball_X;
			//	ball_X = ChSharedPtr<ChBody>(my_system->Get_bodylist()->at(0 + middleBall));
			//	my_system->Get_bodylist()->at(0 + middleBall)->AddRef();

			//	fprintf(PosVelWrite2, "%f,%f,%f,%f,%f,%f,\n",
			//		ball_X->GetPos().x, ball_X->GetPos().y, ball_X->GetPos().z,
			//		ball_X->GetPos_dt().x, ball_X->GetPos_dt().y, ball_X->GetPos_dt().z);
			//}

			//fclose(PosVelWrite2);

			// _________________________________________________________________________
			// fifth row information

			FILE * PosVelWrite2 = fopen(particlePos2_file.c_str(), "a");
			{
				ChSharedPtr<ChBody> ball_5;
				ball_5 = ChSharedPtr<ChBody>(my_system->Get_bodylist()->at(0 + fifthBall));
				my_system->Get_bodylist()->at(0 + fifthBall)->AddRef();

				fprintf(PosVelWrite2, "%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,\n",
				//fprintf(PosVelWrite2, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,\n",
					ball_5->GetPos().x, ball_5->GetPos().y, ball_5->GetPos().z,
					ball_5->GetPos_dt().x, ball_5->GetPos_dt().y, ball_5->GetPos_dt().z);
				printf( "%.10f\t", ball_5->GetPos().y);
			}

			fclose(PosVelWrite2);

			// _________________________________________________________________________
			// first row information

			FILE * PosVelWrite3 = fopen(particlePos3_file.c_str(), "a");
			{
				ChSharedPtr<ChBody> ball_1;
				ball_1 = ChSharedPtr<ChBody>(my_system->Get_bodylist()->at(0 + secondBall));
				my_system->Get_bodylist()->at(0 + secondBall)->AddRef();

				//fprintf(PosVelWrite2, "%.15f,%.15f,%.15f,%.15f,%.15f,%.15f,\n",
				fprintf(PosVelWrite3, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,\n",
					ball_1->GetPos().x, ball_1->GetPos().y, ball_1->GetPos().z,
					ball_1->GetPos_dt().x, ball_1->GetPos_dt().y, ball_1->GetPos_dt().z);
				printf("%.10f\n", ball_1->GetPos().y);
			}

			fclose(PosVelWrite3);

			// _________________________________________________________________________
			//// all balls position information
			//{
			//	FILE * AllWrite2 = fopen(allPos_file.c_str(), "w");
			//	for (int ii = 0; ii < 22; ii++){
			//		ChSharedPtr<ChBody> ball_ii;
			//		ball_ii = ChSharedPtr<ChBody>(my_system->Get_bodylist()->at(ii));
			//		my_system->Get_bodylist()->at(ii)->AddRef();

			//		fprintf(AllWrite2, "%d,%.5e,%.5e,%.5e,%.5e,%.5ef,%.5e,\n",
			//			//fprintf(PosVelWrite2, "%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,\n",
			//			ii,
			//			ball_ii->GetPos().x, ball_ii->GetPos().y, ball_ii->GetPos().z,
			//			ball_ii->GetPos_dt().x, ball_ii->GetPos_dt().y, ball_ii->GetPos_dt().z);
			//		//printf("%.10f\n", ball_5->GetPos().y);
			//	}
			//}

			//fclose(PosVelWrite2);

			//// _________________________________________________________________________
			//// last row information

			//FILE * PosVelWrite3 = fopen(particlePos3_file.c_str(), "a");
			//{
			//	ChSharedPtr<ChBody> ball_X;
			//	ball_X = ChSharedPtr<ChBody>(my_system->Get_bodylist()->at(0 + lastBall - 1));
			//	my_system->Get_bodylist()->at(0 + lastBall - 1)->AddRef();

			//	fprintf(PosVelWrite3, "%f,%f,%f,%f,%f,%f,\n",
			//		ball_X->GetPos().x, ball_X->GetPos().y, ball_X->GetPos().z,
			//		ball_X->GetPos_dt().x, ball_X->GetPos_dt().y, ball_X->GetPos_dt().z);
			//}

			//fclose(PosVelWrite3);
			
			data_out_frame++;
		} // end of if - writing data to file

		// ________________________________________________________________
		// Produce POV-Ray data
		if (write_povray_data && my_system->GetChTime() >= visual_out_frame * visual_out_step) {
			char filename[100];
			sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), visual_out_frame + 1);
			utils::WriteShapesPovray(my_system, filename, false);

			visual_out_frame++;
		}
	}// end of while

	cout << "Write checkpoint data to " << flush;
	utils::WriteCheckpoint(my_system, checkPoint_file);

	return 0;
	
}