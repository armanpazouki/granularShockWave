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
// Authors: Michal Kwarta, Arman Pazouki
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

#include "utils/ChUtilsCreators.h"
#include "utils/ChUtilsGenerators.h"
#include "utils/ChUtilsInputOutput.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/lcp/ChLcpSystemDescriptorParallel.h"

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::flush;
using std::endl;

// =============================================================
// Above headers and namespaces

#define VELOCITY_EXCITATION

#undef CHRONO_OPENGL

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

double geomScale = 1e6;  // 1;

// Parameters you set up to have different simulations
double overlap = 1e-3;            // 1e-9 * geomScale; //
double time_step = 1e-13;         // 1.0e-13;
double velocity = 5 * geomScale;  // excitation velocity || m/s

double simulation_time = 5e-8;
bool write_povray_data = true;
double data_out_step = 50 * time_step; // 1e-8;       // time interval between data outputs
double visual_out_step = 50 * time_step; // 1e-7;     // time interval between PovRay outputs
double energy_out_step = 1000 * time_step;  // interval between output
											// conderning potential and kinetoc
											// energy

enum RestrainModel {
	constraint, force
};
RestrainModel resModel = constraint;

const std::string out_dir1 = "Seattle_3";
const std::string out_dir = out_dir1 + "/2dFreeSurface_101x100_vel=100_noGrav";

const std::string pov_dir = out_dir + "/POVRAY";

const std::string pov_vel_dir = out_dir + "/POVRAY_velocities";
const std::string energy_dir = out_dir + "/Energy";
const std::string particles_dir = out_dir + "/Particles";

const std::string time_file = out_dir + "/time.dat";

const std::string checkPoint_file = out_dir + "/checkPoint.dat";

ChSharedPtr<ChBody> ground;
void OutputData(ChSystemParallelDEM* sys, int out_frame, double time) {
	char filename[100];
	sprintf(filename, "%s/data_%04d.dat", pov_dir.c_str(), out_frame);
	utils::WriteShapesPovray(sys, filename, false);
}

double gravity = 0.0;                              //9.81 * geomScale;  // m/s/s

double ballRad = 0.5;  // 0.5e-6 * geomScale; // micrometer
double ballDiam = 2.0 * ballRad;

float mu = 0.18f;
float ballDensity = 2000.0f / geomScale / geomScale / geomScale; // 1.0f; // kg/m/m/m
float ballMass = ballDensity * 4.0f / 3.0f * CH_C_PI * ballRad * ballRad
		* ballRad;

float E = 73e9 / geomScale;  // 1e9 / geomScale;//1e9;// / Young's modulus. Pa =
							 // N/m/m = kg/s/s/m || to watch out for crazy
							 // vibrations
float nu = 0.17f;            // 0.32f;

ChVector<> inertia = (2.0 / 5.0) * ballMass * ballRad * ballRad
		* ChVector<>(1, 1, 1);

double xLay = 101;
double zLay = 100;
double yLay = 1;

double pushedBall = (xLay - 1) / 2;

// xOverlap - for 2D and 3D only
// yOverlap above - different values for ONE_D_ARRANGEMENT and different for
// HEXAGONAL_CLOSE because:
// in ONE_D_ARRANGEMENT - spheres touch each other in the upright way
// in HEXAGONAL_CLOSE - they touch each other in the "triangle" way (equilateral
// triangle to be clear)
// zOverlap - for 3D only

//	ChVector<> CameraLocation;
//	ChVector<> CameraLookAt;

// add some overlap to here
double xDim = (ballDiam - overlap) * xLay + 0.5 * (ballDiam - overlap);
double yDim = ballDiam;
double zDim = (zLay - 1) * (ballDiam - overlap) * sqrt(3) / 2;

ChVector<> CameraLocation = ChVector<>(xDim / 2.0, 10 * yDim, 0);
ChVector<> CameraLookAt = ChVector<>(xDim / 2.0, 0.01, 0.01);

int threads = 2;

bool thread_tuning = false;

// double time_step = 1.0e-9; // is was set up above
double tolerance = 1.0;  // 0.001;   // Arman, check this
int max_iteration_bilateral = 50;

bool clamp_bilaterals = false;
double bilateral_clamp_speed = 0.1;

double timeStartPushing = 0; // 50*time_step;//3.0e-6; // start time for pushing the wall

int data_out_frame = 0;
int visual_out_frame = 0;
int energy_out_frame = 0;

//#ifndef ONE_D_ARRANGEMENT
//	void AddWall(ChSystemParallelDEM* sys, ChSharedPtr<ChBody> wall,
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
//			ChSharedPtr<ChMaterialSurfaceDEM> mat_ext;
//			mat_ext = ChSharedPtr<ChMaterialSurfaceDEM>(new
// ChMaterialSurfaceDEM);
//			mat_ext->SetYoungModulus(Y_ext);
//			mat_ext->SetPoissonRatio(nu_ext);
//			mat_ext->SetRestitution(COR_ext);
//			mat_ext->SetFriction(mu_ext);
//
//		wall->SetMaterialSurface(mat_ext);
//
//		wall->GetCollisionModel()->ClearModel();
//
//		utils::AddBoxGeometry(wall.get_ptr(), dim, ChVector<>(0, 0, 0),
// ChQuaternion<double>(1, 0, 0, 0), visualization);
//		wall->GetCollisionModel()->SetFamily(1);
//		wall->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
//
//		wall->GetCollisionModel()->BuildModel();
//		sys->AddBody(wall);
//	}
//#endif

// =============================================================================
void SetArgumentsForMbdFromInput(int argc, char* argv[]) {
	if (argc > 1) {
		const char* text = argv[1];
		double vel_m_s = atof(text);

		velocity = vel_m_s * geomScale;
	}
	if (argc > 2) {
		const char* text = argv[2];
		xLay = atoi(text);
	}
	if (argc > 3) {
		const char* text = argv[3];
		yLay = atoi(text);
	}
	if (argc > 4) {
		const char* text = argv[4];
		threads = atoi(text);
	}
}
// =============================================================================

void InitializeSystem(ChSystemParallelDEM * my_system) {
	int max_threads = my_system->GetParallelThreadNumber();
	if (threads > max_threads)
		threads = max_threads;
	my_system->SetParallelThreadNumber(threads);
	omp_set_num_threads(threads);

	my_system->GetSettings()->max_threads = threads;
	my_system->GetSettings()->perform_thread_tuning = thread_tuning;

	my_system->GetSettings()->solver.use_full_inertia_tensor = false;
	my_system->GetSettings()->solver.tolerance = tolerance;
	my_system->GetSettings()->solver.max_iteration_bilateral =
			max_iteration_bilateral;
	my_system->GetSettings()->solver.clamp_bilaterals = clamp_bilaterals;
	my_system->GetSettings()->solver.bilateral_clamp_speed =
			bilateral_clamp_speed;

	my_system->GetSettings()->solver.tangential_displ_mode =
			ChSystemDEM::TangentialDisplacementModel::MultiStep;

	my_system->GetSettings()->collision.bins_per_axis = I3(20, 1, 20);
	my_system->GetSettings()->collision.narrowphase_algorithm =
			NARROWPHASE_HYBRID_MPR;

	my_system->GetSettings()->solver.contact_force_model =
			ChSystemDEM::ContactForceModel::Hertz;
	my_system->GetSettings()->solver.adhesion_force_model =
			ChSystemDEM::AdhesionForceModel::DMT;
	my_system->Set_G_acc(ChVector<>(0, 0, -gravity));

}
// =============================================================================

int CreateOutputFolders() {
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
	return 0;
}
// =============================================================================
void BringToPlane(ChSystemParallelDEM * my_system, chrono::int2 bRange) {
	// projections + information about violation
	// .. and some checking of plane violation
	double planeViolation = 0;
	double planeViolation2 = 0;
	//
	double planeViolationAfter = 0;
	double planeViolationAfter2 = 0;

	int numParticles = bRange.y - bRange.x;
	for (int i = 0; i < numParticles; i++) { // Note: it is assumed that particles are written to the system before the ground
		ChSharedPtr<ChBody> ball_X = my_system->Get_bodylist()->at(
				bRange.x + i);

		planeViolation += ball_X->GetPos().y * ball_X->GetPos().y;
		planeViolation2 += fabs(ball_X->GetPos().y);

		ball_X->SetPos(
				ChVector<double>(ball_X->GetPos().x, 0.0, ball_X->GetPos().z));
		ball_X->SetRot(
				ChQuaternion<double>(ball_X->GetRot().e0, 0.0,
						ball_X->GetRot().e2, 0.0));
		ball_X->SetPos_dt(
				ChVector<double>(ball_X->GetPos_dt().x, 0.0,
						ball_X->GetPos_dt().z));
		ball_X->SetRot_dt(
				ChQuaternion<double>(
						-0.5 * ball_X->GetRot().e2 * ball_X->GetWvel_loc().y,
						0.0,
						0.5 * ball_X->GetRot().e2 * ball_X->GetWvel_loc().y,
						0.0));

		planeViolationAfter += ball_X->GetPos().y * ball_X->GetPos().y;
		planeViolationAfter2 += fabs(ball_X->GetPos().y);
	}

	char planeFile[100];
	sprintf(planeFile, "%s/planeViolation.dat", out_dir.c_str());

	FILE* planeViol = fopen(planeFile, "a");
	fprintf(planeViol, "%.14f,%.14f,%.14f,%.14f,\n", planeViolation,
			sqrt(planeViolation2 / numParticles), planeViolationAfter,
			sqrt(planeViolationAfter2 / numParticles));
	fclose(planeViol);
}
// =============================================================================
// Function to calculate system's energy
double TotalEnergy(ChSystemParallelDEM * my_system, chrono::int2 bRange) {
	double kEpE_Energy = 0;

	// -----------------------------------
	// ---- kinetic and potential energy
	// -----------------------------------

	int numParticles = bRange.y - bRange.x;
	for (int i = 0; i < numParticles; i++) {
		ChSharedPtr<ChBody> body = my_system->Get_bodylist()->at(bRange.x + i);
		ChVector<> rotEnergyComps = 0.5 * body->GetInertiaXX()
				* body->GetWvel_loc() * body->GetWvel_loc();
		double kineticEnergy = 0.5 * body->GetMass()
				* body->GetPos_dt().Length2() + rotEnergyComps.x
				+ rotEnergyComps.y + rotEnergyComps.z;
		ChVector<> potVec = body->GetMass()
				* (my_system->Get_G_acc() * body->GetPos());
		double potentialEnergy = potVec.x + potVec.y + potVec.z;
		kEpE_Energy += kineticEnergy + potentialEnergy;
	}

	// -----------------------------------
	// ---- contact energy
	// -----------------------------------

	double contactEnergy = 0;
	int numContacts =
			my_system->data_manager->host_data.dpth_rigid_rigid.size();
	for (int i = 0; i < numContacts; i++) {
		real depth = -1
				* my_system->data_manager->host_data.dpth_rigid_rigid[i]; // depth is a negative number when the two objects are in contact.
		if (depth < 0) {
			continue;
		}
		chrono::int2 ids =
				my_system->data_manager->host_data.bids_rigid_rigid[i];
		ChSharedPtr<ChBody> body1 = my_system->Get_bodylist()->at(ids.x);
		ChSharedPtr<ChBody> body2 = my_system->Get_bodylist()->at(ids.y);

		real2* elastic_moduli =
				my_system->data_manager->host_data.elastic_moduli.data();
		real Y1 = elastic_moduli[ids.x].x;
		real Y2 = elastic_moduli[ids.y].x;
		real nu1 = elastic_moduli[ids.x].y;
		real nu2 = elastic_moduli[ids.y].y;
		real inv_E = (1 - nu1 * nu1) / Y1 + (1 - nu2 * nu2) / Y2;
		real inv_G = 2 * (2 + nu1) * (1 - nu1) / Y1
				+ 2 * (2 + nu2) * (1 - nu2) / Y2;

		real E_eff = 1 / inv_E;
		real r_eff = my_system->data_manager->host_data.erad_rigid_rigid[i];

		double HertzEnergy = 8.0 / 15.0 * E_eff * chrono::Sqrt(r_eff)
				* pow(depth, 2.5);

		real* adhesionMultDMT =
				my_system->data_manager->host_data.adhesionMultDMT_data.data();
		real adhesionMultDMT_eff = chrono::Min(adhesionMultDMT[ids.x],
				adhesionMultDMT[ids.x]);

		double DMT_Energy = adhesionMultDMT_eff * chrono::Sqrt(r_eff) * depth;

		contactEnergy += (HertzEnergy - DMT_Energy);
	}
	return kEpE_Energy + contactEnergy;
}
// =============================================================================

void CreateGround(ChSystemParallelDEM * my_system) {
	ChSharedPtr<ChMaterialSurfaceDEM> groundMat = ChSharedPtr<
			ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
	ground = ChSharedPtr<ChBody>(new ChBody(new ChCollisionModelParallel));
	ground->SetMaterialSurface(groundMat);
	ground->SetIdentifier(0);
	ground->SetMass(10 * ballMass);
	ground->SetPos(ChVector<>(0, -3 * ballRad, 0)); // spheres are generated at y = 0
	ground->SetCollide(false);
	ground->SetBodyFixed(true);
	ground->GetCollisionModel()->ClearModel();
	utils::AddBoxGeometry(ground.get_ptr(),
			ChVector<>(10 * ballRad, 0.1 * ballRad, 10 * ballRad),
			ChVector<>(0, 0, 0), ChQuaternion<>(1, 0, 0, 0), true);	//end wall
	ground->GetCollisionModel()->BuildModel();
	my_system->AddBody(ground);
}
// =============================================================================

void AddBall(ChSystemParallelDEM* sys, double x, double y, double z, int ballId,
		bool ifFixed) {
	ChSharedPtr<ChMaterialSurfaceDEM> ballMat;
	ballMat = ChSharedPtr<ChMaterialSurfaceDEM>(new ChMaterialSurfaceDEM);
	ballMat->SetYoungModulus(E);
	ballMat->SetPoissonRatio(nu);
	ballMat->SetFriction(mu);
	ballMat->SetRestitution(0);  // needed in DMT model
	ballMat->SetAdhesionMultDMT(4.0 / 3.0 * E * pow(overlap, 1.5));

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
	utils::AddSphereGeometry(ball.get_ptr(), ballRad, ChVector<double>(0, 0, 0),
			ChQuaternion<double>(1, 0, 0, 0), true);
	ball->GetCollisionModel()->SetFamily(1);
	ball->GetCollisionModel()->SetDefaultSuggestedEnvelope(0.00);
	ball->GetCollisionModel()->SetDefaultSuggestedMargin(0.00 * ballDiam);

	ball->GetCollisionModel()->BuildModel();

	sys->AddBody(ball);
}

// =============================================================================

chrono::int2 AddBalls(ChSystemParallelDEM* sys) {
	chrono::int2 range;
	range.x = sys->Get_bodylist()->size();
	int ballId = 0;

	bool ifFixed = false; // before figuring out periodic boundaries let spheres
						  // at the boundaries be fixed

	double x = 0., z = 0.;
	double y = 0.;
	double dimention = ballDiam - overlap;// dimention = length of the triangle's side (it
										  // repeats all the time in the code below - that is
										  // why it is no sense in computing it all the time)

	for (int ii = 0; ii < zLay; ii++) {
		z = ii * (-dimention) * sqrt(3) / 2;  // height of the triangle

		for (int jj = 0; jj < xLay; jj++) {
			// ========= if Fixed ??? =================
			ifFixed = false;

			if (ii == zLay - 1)  // last row should be fixed
				ifFixed = true;

			if ((jj == 0) || (jj == xLay - 1)) // first and last column should be fixed
				ifFixed = true;
			// =========================================

			if (ii % 2 == 0) {
				x = jj * dimention;
				ballId++;
				AddBall(sys, x, y, z, ballId, ifFixed);	// note that ballId starts from 1. Reserve 0 for the ground
			} else {
				x = jj * dimention + dimention * 0.5;
				ballId++;
				AddBall(sys, x, y, z, ballId, ifFixed);	// note that ballId starts from 1. Reserve 0 for the ground
			}
		}
	}
	range.y = sys->Get_bodylist()->size();
	return range;
}
// =============================================================================
// Function to calculate system's energy
void RestrainParticles(ChSystemParallelDEM * my_system, chrono::int2 bRange) {
	int numParticles = bRange.y - bRange.x;
	for (int i = 0; i < numParticles; i++) {
		ChSharedPtr<ChBody> body = my_system->Get_bodylist()->at(bRange.x + i);

		// Constrain spindles in a horizontal plane (based on current post location)
		ChSharedPtr<ChLinkLockPointPlane> part_on_YPlane = ChSharedPtr<
				ChLinkLockPointPlane>(new ChLinkLockPointPlane());
		part_on_YPlane->SetNameString("particle_plane_constraint");
		part_on_YPlane->Initialize(body, ground,
				ChCoordsys<>(body->GetPos(), Q_from_AngAxis(CH_C_PI / 2, VECT_X))); // By default, the constraint is in the z_direction. Therefore, it needs to rotate to the y-direction.
		my_system->AddLink(part_on_YPlane);
	}
}
// =============================================================================
// bRange denotes the range of particles in ChSystem body list; i.e. start and end indices
void PrintResults(ChSystemParallelDEM * my_system, chrono::int2 bRange) {
	static int count = -1;
	count++;
	if (my_system->GetChTime() >= data_out_frame * data_out_step) {
		// _________________________________________________________________________
		// time information
		FILE* timeWrite = fopen(time_file.c_str(), "a");
		fprintf(timeWrite, "%.13f,\n", my_system->GetChTime());
		fclose(timeWrite);

		const int howManyParticles = 6;
		int particleNumber[howManyParticles];
		particleNumber[0] = 0;// the middle one (0th) in the first row (that is why
							  // there is an odd amount of spheres in the first row)
		particleNumber[1] = 10;
		particleNumber[2] = 20;
		particleNumber[3] = 30;
		particleNumber[4] = 40;
		particleNumber[5] = 50;

		// recalculating particles' numbers for their real ID
		const int numOfParticles = 2 * howManyParticles - 1;// -1 because particle 0 is the same particle
															// in row 0... we don't want to have same
															// output twice
		int particleID[numOfParticles];
		for (int ii = 0; ii < howManyParticles; ii++) {
			particleID[ii] = (xLay - 1) / 2 + particleNumber[ii];
		}
		for (int ii = 0; ii < howManyParticles - 1; ii++) { // note that number of iterations in the loop changed
			particleID[ii + howManyParticles] = (particleNumber[ii + 1] - 1)
					* xLay + (xLay - 1) / 2;
		}

		// _________________________________________________________________________
		// information about certain particles' position and velocities

		for (int jj = 0; jj < numOfParticles; jj++) {
			char rowFileName[100];
			sprintf(rowFileName, "%s/particle%04d.dat", particles_dir.c_str(),
					particleID[jj]);
			FILE* rowFile = fopen(rowFileName, "a");

			ChSharedPtr<ChBody> ball_X = my_system->Get_bodylist()->at(
					bRange.x + particleID[jj]);

			fprintf(rowFile, "%.13f,%.13f,%.13f,%.13f,%.13f,%.13f,%.13f,\n",
					my_system->GetChTime(), ball_X->GetPos().x,
					ball_X->GetPos().z, ball_X->GetRotAngle(),
					ball_X->GetPos_dt().x, ball_X->GetPos_dt().z,
					ball_X->GetWvel_par().y);

			fclose(rowFile);
		}
		data_out_frame++;
	}  // end of if - writing data to file

// ________________________________________________________________
// print energy
	const std::string fileEnergy = std::string("Energy.csv");
	std::ofstream outEnergy;
	if (count == 0) {
		outEnergy.open(fileEnergy.c_str());
	} else {
		outEnergy.open(fileEnergy.c_str(), std::ios::app);
	}
	outEnergy << my_system->GetChTime() << ", "
			<< TotalEnergy(my_system, bRange) << std::endl;
// ________________________________________________________________
// Produce POV-Ray data
	if (write_povray_data
			&& my_system->GetChTime() >= visual_out_frame * visual_out_step) {
		// =====================================================================
		// =====================================================================

		char filename[100];
		sprintf(filename, "%s/data_%04d.dat", pov_dir.c_str(),
				visual_out_frame + 1);
		utils::WriteShapesPovray(my_system, filename, false);

		// =====================================================================
		// =====================================================================

		// ... and velocity magnitude of spheres
		char AllVel_file[100];
		sprintf(AllVel_file, "%s/data_%04d.dat", pov_vel_dir.c_str(),
				visual_out_frame + 1);
		FILE* AllVel_a = fopen(AllVel_file, "w");

		double maxVelocity = 0;

		for (int i = bRange.x; i < bRange.y; i++) {
			ChSharedPtr<ChBody> ball_X = ChSharedPtr<ChBody>(
					my_system->Get_bodylist()->at(i));

			double velo = sqrt(
					ball_X->GetPos_dt().x * ball_X->GetPos_dt().x
							+ ball_X->GetPos_dt().y * ball_X->GetPos_dt().y
							+ ball_X->GetPos_dt().z * ball_X->GetPos_dt().z);

			fprintf(AllVel_a, "%f,\n", velo);

			if (maxVelocity < velo) {
				maxVelocity = velo;
			}
		}

		fclose(AllVel_a);

		// =====================================================================
		// =====================================================================

		// maximal velocity
		char MaxVelocity_file[100];
		sprintf(MaxVelocity_file, "%s/maxVelocity_%04d.dat",
				pov_vel_dir.c_str(), visual_out_frame + 1);
		FILE* MaxVelocity_a = fopen(MaxVelocity_file, "w");
		fprintf(MaxVelocity_a, "%f,\n", maxVelocity);
		fclose(MaxVelocity_a);

		// =====================================================================
		// =====================================================================

		// ... colors - just to have sth to color them with
		char AllVelColor_file[100];
		sprintf(AllVelColor_file, "%s/dataColor_%04d.dat", pov_vel_dir.c_str(),
				visual_out_frame + 1);
		FILE* AllVelColor_a = fopen(AllVelColor_file, "w");

		ChSharedPtr<ChBody> ball_X = ChSharedPtr<ChBody>(
				my_system->Get_bodylist()->at(bRange.y - 1));
		fprintf(AllVelColor_a, "%.9f,\n",
				ball_X->GetPos_dt().Length() / velocity); // Arman: I don't know what is this, just kept it here for now.

		for (int i = bRange.x; i < bRange.y; i++) {
			ChSharedPtr<ChBody> ball_X = ChSharedPtr<ChBody>(
					my_system->Get_bodylist()->at(i));
			double velo = ball_X->GetPos_dt().Length();
			fprintf(AllVelColor_a, "%.9f,\n", velo / velocity);
		}

		fclose(AllVelColor_a);
		visual_out_frame++;
	}
}
// =============================================================================

int main(int argc, char* argv[]) {
	// ================================================
	// ================================================
	SetArgumentsForMbdFromInput(argc, argv);
	ChSystemParallelDEM* my_system = new ChSystemParallelDEM();
	InitializeSystem(my_system);

	const int howManyParticles = 6;
	int particleNumber[howManyParticles];
	particleNumber[0] = 0;// the middle one (0th) in the first row (that is why
						  // there is an odd amount of spheres in the first row)
	particleNumber[1] = 10;
	particleNumber[2] = 20;
	particleNumber[3] = 30;
	particleNumber[4] = 40;
	particleNumber[5] = 50;

	// recalculating particles' numbers for their real ID
	const int numOfParticles = 2 * howManyParticles - 1;// -1 because particle 0 is the same particle
														// in row 0... we don't want to have same
														// output twice
	int particleID[numOfParticles];
	for (int ii = 0; ii < howManyParticles; ii++) {
		particleID[ii] = (xLay - 1) / 2 + particleNumber[ii];
	}
	for (int ii = 0; ii < howManyParticles - 1; ii++) { // note that number of iterations in the loop changed
		particleID[ii + howManyParticles] = (particleNumber[ii + 1] - 1) * xLay
				+ (xLay - 1) / 2;
	}

	// ================================================
	// ================================================

	// cout << "ballDensity" << ballDensity;

	int isGood = CreateOutputFolders();
	if (isGood == 1) {
		return 1;
	}

	cout << "Create DEM-P system" << endl;

#ifdef CHRONO_OPENGL
	opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
	const std::string title = "Wave propagation in granular media. DEM-P.";
	gl_window.Initialize(800, 600, title.c_str(), my_system);

	gl_window.SetCamera(CameraLocation, ChVector<>(0.2, 0.3, 0.4),
			ChVector<>(0, 0, 1), ballDiam, ballDiam);
//	gl_window.viewer->render_camera.camera_scale = 2.0/(1000.0);
//	gl_window.viewer->render_camera.near_clip = .001;
// gl_window.SetRenderMode(opengl::WIREFRAME);

#endif

	chrono::int2 bRange = AddBalls(my_system);
	CreateGround(my_system); // preferably, call CreateGround after AddBalls due to indexing issues of the balls in body_list; although I must have figured it already.
	if (resModel == constraint) {
		RestrainParticles(my_system, bRange);
	}

	bool pushed = false;

	while (my_system->GetChTime() <= simulation_time) {
#ifdef VELOCITY_EXCITATION
		//// push the first sphere after the certain time
		if (pushed == false && my_system->GetChTime() >= timeStartPushing) {
			pushed = true;
			ChSharedPtr<ChBody> ball_X = ChSharedPtr<ChBody>(
					my_system->Get_bodylist()->at(bRange.x + pushedBall));
			ball_X->SetPos_dt(
					ball_X->GetPos_dt() + ChVector<>(0, 0, -velocity));
		}
#else
// Force excitation
#endif

#ifdef CHRONO_OPENGL
		if (gl_window.Active()) {
			gl_window.DoStepDynamics(time_step);
			gl_window.Render();
		} else {
			break;
		}
#else
		my_system->DoStepDynamics(time_step);
#endif

		if (resModel == force) {
			BringToPlane(my_system, bRange);
		}
		PrintResults(my_system, bRange);
		if (my_system->GetChTime() >= energy_out_frame * energy_out_step) {

		}

	}  // end of while

	cout << "Write checkpoint data to " << flush;
	utils::WriteCheckpoint(my_system, checkPoint_file);

	return 0;
}
