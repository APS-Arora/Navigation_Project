#include "c_Integration.h"

#include <iostream>
#include <string>
#include <random>
#include <chrono>
#include <time.h>
#include <cmath>
#include <fstream>
#include <vector>
#include <Eigen\Dense>
#include "c_main.h"


CInt::CInt()
{
}


CInt::~CInt()
{
}

void CInt::m_predict(INS::INS_States INS_Estimate, CMain::InsOutput IMU_Output)
{
	Matrix3d Omega_earth = m_Skew(Vector3d::UnitZ()*RtEarthRotn_Const);

	Vector3d corr_sf = IMU_Output.sf_actual - INS_Estimate.accel_bias;
	Vector3d corr_angvel = IMU_Output.ang_actual - INS_Estimate.gyro_bias;

	// Propagating State
	m_ErrorState.clock_offset += m_ErrorState.clock_drift*dt;

	// Calculating State Transition Matrix
	Matrix3d Fe21 = -m_Skew(INS_Estimate.Cb_e*corr_sf);
	Matrix3d Fe23;
	Matrix<double, 17, 17> Phi = Matrix<double,17,17>::Identity();

}

Matrix3d CInt::m_Skew(Vector3d a){
	Matrix3d A;
	A << 0, -a(2), a(1),
		a(2), 0, -a(0),
		-a(1), a(0), 0;
	return A;
}