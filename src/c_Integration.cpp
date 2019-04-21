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

void CInt::InitErrorCov(){
	double mg_ms2 = 0.00000980665, deg_to_rad = 0.01745329252;
	ErrorCov.topLeftCorner(3, 3) = Matrix3d::Identity()*pow(0.0001, 2);
	ErrorCov.block(3, 3, 3, 3) = Matrix3d::Identity()*pow(0.1, 2);
	ErrorCov.block(3, 3, 6, 6) = Matrix3d::Identity()*pow(10, 2);
	ErrorCov.block(3, 3, 9, 9) = Matrix3d::Identity()*pow(1000 * mg_ms2, 2);
	ErrorCov.block(3, 3, 12, 12) = Matrix3d::Identity()*pow(10 * deg_to_rad / 3600, 2);
	ErrorCov(15, 15) = pow(10, 2);
	ErrorCov(16, 16) = pow(0.1, 2);
}

void CInt::m_predict(CMain::INS_States INS_Estimate, CMain::InsOutput IMU_Output)
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

void CInt::LsPosVel(CMain::SatData *mp_GpsSatData, CMain::INS_States& INS_Init){

	//Speed of Light
	double c = 299792458, omega_ie = 0.00007292115;
	Vector4d x_pred = { 0, 0, 0, 0 }, x_est;
	Matrix<double, Dynamic, 4> H_mat = Matrix<double, Dynamic, 4>::Zero(no_sat, 4);
	Matrix<double, Dynamic, 3> gnss_pos = Matrix<double, Dynamic, 3>::Zero(no_sat, 3),
		gnss_vel = Matrix<double, Dynamic, 3>::Zero(no_sat, 3);
	VectorXd pred_meas = VectorXd::Zero(no_sat), gnss_range = VectorXd::Zero(no_sat), gnss_range_rate = VectorXd::Zero(no_sat);
	Vector3d delta_r, est_pos, omega, u_as_e;
	Matrix3d Ce_i, Omega_ie;
	double approx_range, range, range_rate;
	int cnvg = 1;

	while (cnvg > 0.0001){
		for (int j = 0; j < no_sat; j++){
			//Check if satellite is visible
			if (!mp_GpsSatData[j].visible){
				continue;
			}
			gnss_pos(j, 0) = mp_GpsSatData[j].x_cord; gnss_pos(j, 1) = mp_GpsSatData[j].y_cord, gnss_pos(j, 2) = mp_GpsSatData[j].z_cord;
			gnss_range(j) = mp_GpsSatData[j].range;
			//Predict approx range
			delta_r = gnss_pos.col(j) - x_pred.head(3);
			approx_range = delta_r.norm();

			//Frame rotation during signal transit time (8.36)
			Ce_i << 1, omega_ie * approx_range / c, 0,
				-omega_ie * approx_range / c, 1, 0,
				0, 0, 1;

			//Predict pseudorange (9.144)
			delta_r = Ce_i*gnss_pos.col(j) - x_pred.head(3);
			range = delta_r.norm();
			pred_meas[j] = range + x_pred(4);

			//Predict line of sight and deploy in measurement matrix, (9.144)
			H_mat.block<1, 3>(j, 0) = -delta_r.transpose() / range;
			H_mat(j, 4) = 1;
		}
		//Unweighted least - squares solution
		x_est = x_pred + (H_mat.transpose()*H_mat).inverse()*H_mat.transpose()*(gnss_range - pred_meas);

		//Test convergence
		cnvg = (x_est - x_pred).norm();

		//Set predictions to estimates for next iteration
		x_pred = x_est;
	}
	INS_Init.position = x_est.head(3);
	m_ErrorState.clock_offset = x_est(4);

	//VELOCITY AND CLOCK DRIFT
	omega << 0, 0, omega_ie;
	Omega_ie = m_Skew(omega);

	//Resetting values
	x_pred = { 0, 0, 0, 0 }, cnvg = 1;

	while (cnvg > 0.0001){
		for (int j = 0; j < no_sat; j++){
			//Check if satellite is visible
			if (!mp_GpsSatData[j].visible){
				continue;
			}
			gnss_vel(j, 0) = mp_GpsSatData[j].x_vel; gnss_vel(j, 0) = mp_GpsSatData[j].y_vel; gnss_vel(j, 0) = mp_GpsSatData[j].z_vel;
			gnss_range_rate(j) = mp_GpsSatData[j].range_rt;
			//Predict approx range
			delta_r = gnss_pos.col(j) - x_pred.head(3);
			approx_range = sqrt(delta_r.transpose()*delta_r);

			//Frame rotation during signal transit time (8.36)
			Ce_i << 1, omega_ie * approx_range / c, 0,
				-omega_ie * approx_range / c, 1, 0,
				0, 0, 1;

			//Predict pseudorange (9.144)
			delta_r = Ce_i*gnss_pos.col(j) - x_pred.head(3);
			range = delta_r.norm();

			//Calculate line of sight using (8.41)
			u_as_e = delta_r / range;

			//Predict range rate
			range_rate = u_as_e.transpose()*(Ce_i*gnss_vel.col(j) + omega_ie*gnss_pos.col(j) - x_pred.head(3) + Omega_ie*INS_Init.position);
			pred_meas(j) = range_rate + x_pred(4);

			//Predict line of sight and deploy in measurement matrix
			H_mat.block<1, 3>(j, 0) = -u_as_e.transpose();
			H_mat(j, 4) = 1;
		}
		//Unweighted least-squares solution, (9.35)/(9.141)
		x_est = x_pred + (H_mat.transpose()*H_mat).inverse()*H_mat.transpose()*(gnss_range_rate - pred_meas);

		//Test convergence
		cnvg = sqrt((x_est - x_pred).transpose()*(x_est - x_pred));

		//Set predictions to estimates for next iteration
		x_pred = x_est;
	}

	//Set outputs to estimates
	INS_Init.velocity = x_est.head(3);
	m_ErrorState.clock_drift = x_est(4);
}

void CInt::m_correct(CMain::INS_States& INS_Estimate, CMain::GNSS_Measurement* GPS_Output)
{
	MatrixXd u_as_e_T = MatrixXd::Zero(no_sat, 3), pred_meas = MatrixXd::Zero(no_sat, 2);
	VectorXd delta_z = VectorXd::Zero(2*no_sat);
	Matrix3d C_e_I, Ce_n=INS_Estimate.Cb_n*INS_Estimate.Cb_e.inverse();
	Matrix3d Omega_earth = m_Skew(Vector3d::UnitZ()*RtEarthRotn_Const);
	Vector3d delta_r, LOS_ned;
	double range, range_rate, azimuth, elevation;
	CMain::AtmSatDelay delays;
	int visible_sats = 0;

	for (int j = 0; j < no_sat; j++)
	{
		if (!GPS_Output[j].visible)
			continue;

		// Predict approx range
		delta_r = GPS_Output[j].SatPos - INS_Estimate.position;
		range = delta_r.norm();

		// Calculate frame rotation during signal transit time
		C_e_I << 1, RtEarthRotn_Const*range / SpeedLight_Const, 0, -RtEarthRotn_Const*range / SpeedLight_Const, 1, 0, 0, 0, 1;

		// Predict range
		delta_r = C_e_I*GPS_Output[j].SatPos - INS_Estimate.position;
		range = delta_r.norm();

		// Predict LOS
		u_as_e_T.row(visible_sats)=delta_r.transpose()/range;

		// Predict Azimuth Elevation
		LOS_ned = Ce_n*delta_r.normalized();
		azimuth = atan2(LOS_ned(1), LOS_ned(0));
		elevation = -asin(LOS_ned(2));

		// Compute Delays
		m_Delay.ComputeDelays(GPS_Output[j].time_data.time_of_week,
			GPS_Output[j].time_data.days_in_year,
			INS_Estimate.latitude,
			INS_Estimate.longitude,
			INS_Estimate.height,
			azimuth, elevation, delays, m_DelayParams);

		// Predict pseudo-range
		pred_meas(visible_sats, 0) = range + delays.tropo_delay_hop + delays.iono_delay_klob - GPS_Output[j].clock_correction + m_ErrorState.clock_offset;
		delta_z(visible_sats) = GPS_Output[j].pseudo_range - pred_meas(visible_sats, 0);

		// Predict range_rate and pseudo-range_rate
		range_rate = delta_r.normalized().dot(C_e_I*(GPS_Output[j].SatVel + Omega_earth*GPS_Output[j].SatPos) - (INS_Estimate.velocity + Omega_earth*INS_Estimate.position));
		pred_meas(visible_sats, 1) = range_rate + m_ErrorState.clock_drift;
		delta_z(visible_sats+no_sat) = GPS_Output[j].pseudo_range_rate - pred_meas(visible_sats, 1);
		visible_sats++;
	}

	// Measurement Matrix
	MatrixXd H = MatrixXd::Zero(2 * visible_sats, 17);
	H.block(0, 6, visible_sats, 3) = u_as_e_T.block(0, 0, visible_sats, 3);
	H.block(0, 15, visible_sats, 1) = VectorXd::Ones(visible_sats);
	H.block(visible_sats, 3, visible_sats, 3) = u_as_e_T.block(0, 0, visible_sats, 3);
	H.block(visible_sats, 16, visible_sats, 1) = VectorXd::Ones(visible_sats);

	// Measurement Noise Covariance
	MatrixXd R = (VectorXd(2 * visible_sats) << m_NoiseConfig.pseudo_range_PSD*VectorXd::Ones(visible_sats), m_NoiseConfig.pseudo_range_rate_PSD*VectorXd::Ones(visible_sats)).finished().asDiagonal();

	// Kalman Gain
	MatrixXd K = ErrorCov*H.transpose()*(H*ErrorCov*H.transpose() + R).inverse();

	// Measurement Innovation
	delta_z << delta_z.head(visible_sats), delta_z.segment(no_sat, visible_sats);

	// State Update
	VectorXd delta_x = K*delta_z;
	m_ErrorState.attitude += delta_x.head(3);
	m_ErrorState.velocity += delta_x.segment(3, 3);
	m_ErrorState.position += delta_x.segment(6, 3);
	m_ErrorState.accel_bias += delta_x.segment(9, 3);
	m_ErrorState.gyro_bias += delta_x.segment(12, 3);
	m_ErrorState.clock_offset += delta_x(15);
	m_ErrorState.clock_drift += delta_x(16);

	// State Covariance Update
	ErrorCov = (MatrixXd::Identity(17, 17) - K*H)*ErrorCov;

	// Close Loop Corrections
	INS_Estimate.Cb_e = (Matrix3d::Identity() - m_Skew(m_ErrorState.attitude))*INS_Estimate.Cb_e;
	m_ErrorState.attitude = Vector3d::Zero();
	INS_Estimate.velocity = INS_Estimate.velocity - m_ErrorState.velocity;
	m_ErrorState.velocity = Vector3d::Zero();
	INS_Estimate.position = INS_Estimate.position - m_ErrorState.position;
	m_ErrorState.position = Vector3d::Zero();
	INS_Estimate.accel_bias = INS_Estimate.accel_bias - m_ErrorState.accel_bias;
	m_ErrorState.accel_bias = Vector3d::Zero();
	INS_Estimate.gyro_bias = INS_Estimate.gyro_bias - m_ErrorState.gyro_bias;
	m_ErrorState.gyro_bias = Vector3d::Zero();
}
