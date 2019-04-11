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
	double R_0 = EqtRadiusOfEarth_Const, e = Eccentricity_Const;

	Vector3d corr_sf = IMU_Output.sf_actual - INS_Estimate.accel_bias;
	Vector3d corr_angvel = IMU_Output.ang_actual - INS_Estimate.gyro_bias;

	// Propagating State
	m_ErrorState.clock_offset += m_ErrorState.clock_drift*dt;

	// Calculating State Transition Matrix
	Matrix3d Fe21 = -m_Skew(INS_Estimate.Cb_e*corr_sf);
	double geocentric_radius = R_0 / sqrt(1 - pow(e * sin(INS_Estimate.latitude), 2.)) *sqrt(pow(cos(INS_Estimate.latitude), 2.) + pow(1 - e*e,2.) * pow(sin(INS_Estimate.latitude), 2.));
	Matrix3d Fe23 = 2 * GravityECEF(INS_Estimate.position) / geocentric_radius*INS_Estimate.position.normalized().transpose();
	Matrix<double, 17, 17> Phi;
	Phi << -Omega_earth, MatrixXd::Zero(3, 9), INS_Estimate.Cb_e, MatrixXd::Zero(3, 2), Fe21, -2 * Omega_earth, Fe23, INS_Estimate.Cb_e, Matrix3d::Zero(), MatrixXd::Zero(3, 2), Matrix3d::Zero(), Matrix3d::Identity(), MatrixXd::Zero(3, 17), MatrixXd::Zero(6, 17), MatrixXd::Zero(2, 15), Vector2d(1,0).asDiagonal();
	Phi = Matrix<double, 17, 17>::Identity() + Phi*dt;

	// Calculating Process Noise Covariance
	Matrix<double, 17, 17> Q = (VectorXd(17) << m_NoiseConfig.gyro_PSD, m_NoiseConfig.gyro_PSD, m_NoiseConfig.gyro_PSD, m_NoiseConfig.accel_PSD, m_NoiseConfig.accel_PSD, m_NoiseConfig.accel_PSD, 0, 0, 0, m_NoiseConfig.accel_bias_PSD, m_NoiseConfig.accel_bias_PSD, m_NoiseConfig.accel_bias_PSD, m_NoiseConfig.gyro_bias_PSD, m_NoiseConfig.gyro_bias_PSD, m_NoiseConfig.gyro_bias_PSD, m_NoiseConfig.clock_offset_PSD, m_NoiseConfig.clock_drift_PSD).finished().asDiagonal();

	// Propagating Error Covariance Matrix
	ErrorCov = Phi * (ErrorCov + 0.5 * Q) * Phi.transpose() + 0.5 * Q;
}

Matrix3d CInt::m_Skew(Vector3d a){
	Matrix3d A;
	A << 0, -a(2), a(1),
		a(2), 0, -a(0),
		-a(1), a(0), 0;
	return A;
}

Vector3d CInt::GravityECEF(Vector3d position){
	Vector3d g, vec, gamma;
	double J_2 = 0.001082627;
	double z_scale, mag_r;
	mag_r = sqrt(position.transpose()*position);
	if (!mag_r){
		g << 0, 0, 0;
	}
	else {
		z_scale = 5 * pow(position[2] / mag_r, 2);
		vec << (1 - z_scale)*position[0], (1 - z_scale)*position[1], (3 - z_scale)*position[2];
		gamma = -GMProd_Const / pow(mag_r, 3)*(position + 1.5*J_2*pow(EqtRadiusOfEarth_Const / mag_r, 2)*vec);
		g[0] = gamma[0] + RtEarthRotn_Const*position[0];
		g[1] = gamma[1] + RtEarthRotn_Const*position[1];
		g[2] = gamma[2];
	}
	return g;
}