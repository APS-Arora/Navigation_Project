/********************************************************************************************************
* Description            : Implementation of Tightly Coupled Error State Kalman Filter
* Specific library calls :	1) Using CDelayCalc for delay calculations
* Classes                : COtherSensors
* Assumptions            :  1) Satellite Position and Velocity must be provided to LsPosiVel and Correct. This is to be replaced by ephemeris in future with satellite position, vel calculations to be done inside the class.
							2) LsPosiVel does not correct pseudoranges for tropo and iono delay. Correct does it though.
* Reference              : 1. Chapter 12, Paul D Grooves : Principles of GNSS, Inertial, and Multisensor Integrated Navigation System(2008)
*
* Version History        :
* <1.1><Amanpreetsingh><02/03/2019>
***********************************************************************************************************/
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
#include <unsupported/Eigen/MatrixFunctions>


CInt::CInt(CMain::DelayCalcParam InputDelayParams)
{
	once = true;
	no_sat = 31;

	m_NoiseConfig.accel_PSD = pow(200 * 9.80665e-6,2.0);
	m_NoiseConfig.gyro_PSD = pow((0.02 * EIGEN_PI/180 / 60), 2.0);
	m_NoiseConfig.accel_bias_PSD = 1.0E-7;
	m_NoiseConfig.gyro_bias_PSD = 2.0E-12;
	m_NoiseConfig.clock_drift_PSD = 1;
	m_NoiseConfig.clock_offset_PSD = 1;
	m_NoiseConfig.pseudo_range_PSD = 2.5;
	m_NoiseConfig.pseudo_range_rate_PSD = 0.1;

	m_ErrorState.attitude = Vector3d::Zero();
	m_ErrorState.velocity = Vector3d::Zero();
	m_ErrorState.position = Vector3d::Zero();
	m_ErrorState.accel_bias = Vector3d::Zero();
	m_ErrorState.gyro_bias = Vector3d::Zero();

	dt = 0.1;
	m_DelayParams = InputDelayParams;
}


CInt::~CInt()
{
}

void CInt::InitErrorCov(){
	double mg_ms2 = 0.00000980665, deg_to_rad = 0.01745329252;
	ErrorCov = MatrixXd::Zero(17, 17);
	ErrorCov.topLeftCorner(3, 3) = Matrix3d::Identity()*pow(0.0001, 2);
	ErrorCov.block(3, 3, 3, 3) = Matrix3d::Identity()*pow(0.1, 2);
	ErrorCov.block(6, 6, 3, 3) = Matrix3d::Identity()*pow(10, 2);
	ErrorCov.block(9, 9, 3, 3) = Matrix3d::Identity()*pow(1000 * mg_ms2, 2);
	ErrorCov.block(12, 12, 3, 3) = Matrix3d::Identity()*pow(10 * deg_to_rad / 3600, 2);
	ErrorCov(15, 15) = pow(10, 2);
	ErrorCov(16, 16) = pow(0.1, 2);
}
/********************************************************************************************************
* Function               : Predict
* Description            : Prediction Step for Error State Kalman Filter
* Function Parameter     : None
* Return value           : None
* Specific library calls : None
* Functions called       : None
* Assumptions            : None
* Reference              : None
* Version History        :
* <1.1><Amanpreetsingh><02/03/2019>
***********************************************************************************************************/
void CInt::Predict(CMain::INS_States INS_Estimate, CMain::InsOutput IMU_Output)
{
	Matrix3d Omega_earth = m_Skew(Vector3d::UnitZ()*RtEarthRotn_Const);
	double R_0 = EqtRadiusOfEarth_Const, e = Eccentricity_Const;

	Vector3d corr_sf = IMU_Output.sf_actual - INS_Estimate.accel_bias;
	Vector3d corr_angvel = IMU_Output.ang_actual - INS_Estimate.gyro_bias;

	// Propagating State
	m_ErrorState.clock_offset += m_ErrorState.clock_drift*dt;

	// Calculating State Transition Matrix
	double geocentric_radius = R_0 / sqrt(1 - pow(e * sin(INS_Estimate.latitude), 2.)) *sqrt(pow(cos(INS_Estimate.latitude), 2.) + pow(1 - e*e, 2.) * pow(sin(INS_Estimate.latitude), 2.));
	Matrix<double, 17, 17> Phi = Matrix<double, 17, 17>::Zero();
	Phi.block(0, 0, 3, 3) = -Omega_earth;
	Phi.block(0, 12, 3, 3) = INS_Estimate.Cb_e;
	Phi.block(3, 0, 3, 3) = -m_Skew(INS_Estimate.Cb_e*corr_sf);
	Phi.block(3, 3, 3, 3) = -2 * Omega_earth;
	Phi.block(3, 6, 3, 3) = 2 * GravityECEF(INS_Estimate.position) / geocentric_radius*INS_Estimate.position.normalized().transpose();
	Phi.block(3, 9, 3, 3) = INS_Estimate.Cb_e;
	Phi.block(6, 3, 3, 3) = Matrix3d::Identity();
	Phi(15, 16) = 1;
	/*Phi << -Omega_earth, MatrixXd::Zero(3, 9), INS_Estimate.Cb_e, MatrixXd::Zero(3, 2), Fe21, -2 * Omega_earth, Fe23, INS_Estimate.Cb_e, 
			Matrix3d::Zero(), MatrixXd::Zero(3, 2), Matrix3d::Zero(), Matrix3d::Identity(), MatrixXd::Zero(3, 11), MatrixXd::Zero(6, 17), 
			MatrixXd::Zero(2, 15), Vector2d(1,0).asDiagonal();*/
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
		g[0] = gamma[0] + pow(RtEarthRotn_Const,2.0)*position[0];
		g[1] = gamma[1] + pow(RtEarthRotn_Const,2.0)*position[1];
		g[2] = gamma[2];
	}
	return g;
}

void CInt::LsPosVel(CMain::GNSS_Measurement *mp_GpsSatData, CMain::INS_States& INS_Init){

	//Speed of Light
	double c = 299792458, omega_ie = 0.00007292115, azimuth, elevation;
	Vector4d x_pred = { 0, 0, 0, 0 }, x_est;
	Matrix<double, Dynamic, 4> H_mat = Matrix<double, Dynamic, 4>::Zero(no_sat, 4);
	Matrix<double, 3, Dynamic> gnss_pos = Matrix<double, 3, Dynamic>::Zero(3, no_sat),
		gnss_vel = Matrix<double, 3, Dynamic>::Zero(3, no_sat);
	VectorXd pred_meas = VectorXd::Zero(no_sat), gnss_range = VectorXd::Zero(no_sat), gnss_range_rate = VectorXd::Zero(no_sat);
	Vector3d delta_r, est_pos, omega, u_as_e;//, LOS_ned;
	Matrix3d Ce_i, Omega_ie; //Ce_n = INS_Init.Cb_n*INS_Init.Cb_e.inverse();
	//CMain::AtmSatDelay delays;
	double approx_range, range, range_rate;
	int cnvg = 1;

	while (cnvg > 0.0001){
		for (int j = 0; j < no_sat; j++){
			//Check if satellite is visible
			if (!mp_GpsSatData[j].visible){
				continue;
			}
			gnss_pos.col(j) = mp_GpsSatData[j].SatPos;
			gnss_range(j) = mp_GpsSatData[j].pseudo_range;

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

			/*
			// Predict Azimuth Elevation
			LOS_ned = Ce_n*delta_r.normalized();
			azimuth = atan2(LOS_ned(1), LOS_ned(0));
			elevation = -asin(LOS_ned(2));

			// Compute Delays
			m_Delay.ComputeDelays(mp_GpsSatData[j].time_data.time_of_week,
				mp_GpsSatData[j].time_data.days_in_year,
				INS_Init.latitude,
				INS_Init.longitude,
				INS_Init.height,
				azimuth, elevation, delays, m_DelayParams);
			*/

			pred_meas[j] = range + x_pred(3);

			//pred_meas[j] = range + x_pred(4) + delays.iono_delay_klob + delays.tropo_delay_hop - mp_GpsSatData[j].clock_correction;

			//Predict line of sight and deploy in measurement matrix, (9.144)
			H_mat.block<1, 3>(j, 0) = -delta_r.transpose() / range;
			H_mat(j, 3) = 1;
		}
		//Unweighted least - squares solution
		x_est = x_pred + (H_mat.transpose()*H_mat).inverse()*H_mat.transpose()*(gnss_range - pred_meas);

		//Test convergence
		cnvg = (x_est - x_pred).norm();

		//Set predictions to estimates for next iteration
		x_pred = x_est;
	}
	INS_Init.position = x_est.head(3);
	m_ErrorState.clock_offset = x_est(3);

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
			gnss_vel.col(j) = mp_GpsSatData[j].SatVel;
			gnss_range_rate(j) = mp_GpsSatData[j].pseudo_range_rate;
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
			pred_meas(j) = range_rate + x_pred(3);

			//Predict line of sight and deploy in measurement matrix
			H_mat.block<1, 3>(j, 0) = -u_as_e.transpose();
			H_mat(j, 3) = 1;
		}
		//Unweighted least-squares solution, (9.35)/(9.141)
		x_est = x_pred + (H_mat.transpose()*H_mat).inverse()*H_mat.transpose()*(gnss_range_rate - pred_meas);

		//Test convergence
		cnvg = (x_est - x_pred).norm();

		//Set predictions to estimates for next iteration
		x_pred = x_est;
	}

	//Set outputs to estimates
	INS_Init.velocity = x_est.head(3);
	m_ErrorState.clock_drift = x_est(3);
}

void CInt::run(CMain::INS_States& INS_Estimates, CMain::GNSS_Measurement* GPS_Output, CMain::InsOutput IMU_Output)
{
	if (once)
	{
		LsPosVel(GPS_Output, INS_Estimates);
		Calc_NED_States(INS_Estimates);
		InitAttitude(&IMU_Output, INS_Estimates);
		Calc_NED_States(INS_Estimates);
		InitErrorCov();
		INS_Estimates.accel_bias = INS_Estimates.gyro_bias = Vector3d::Zero();
		once = false;
	}
	INS_Estimate(INS_Estimates, &IMU_Output);
	Predict(INS_Estimates, IMU_Output);
	Correct(INS_Estimates, GPS_Output);
	Calc_NED_States(INS_Estimates);
}
/********************************************************************************************************
* Function               : Correct
* Description            : Correction Step for Error State Kalman Filter
* Function Parameter     : None
* Return value           : None
* Specific library calls : None
* Functions called       : None
* Assumptions            : None
* Reference              : None
* Version History        :
* <1.1><Amanpreetsingh><02/03/2019>
***********************************************************************************************************/
void CInt::Correct(CMain::INS_States& INS_Estimate, CMain::GNSS_Measurement* GPS_Output)
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
	delta_z = (VectorXd(2 * visible_sats) << delta_z.head(visible_sats), delta_z.segment(no_sat, visible_sats)).finished();

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
/********************************************************************************************************
* Function               : INS_Estimate
* Description            : Implementation of Inertial Navigation Equations
* Function Parameter     : None
* Return value           : None
* Specific library calls : None
* Functions called       : Calc_NED_States
* Assumptions            : Implemented in ECEF Reference
* Reference              : None
* Version History        :
* <1.1><Amanpreetsingh><02/03/2019>
***********************************************************************************************************/
void CInt::INS_Estimate(CMain::INS_States& m_INS_States, CMain::InsOutput* m_InsUserOutput)
{
	//
	double alpha_ie = RtEarthRotn_Const*dt;
	Matrix3d C_Earth = m_TransformMatrix(0, 0, alpha_ie);

	//
	Vector3d alpha_ib_b = m_InsUserOutput->ang_actual * dt;
	double mag_alpha = alpha_ib_b.norm();
	Matrix3d Alpha_ib_b = m_Skew(alpha_ib_b);

	//
	Matrix3d C_new_old;
	if (mag_alpha > 1e-8)
	{
		C_new_old = Matrix3d::Identity() + sin(mag_alpha) / mag_alpha*Alpha_ib_b + (1 - cos(mag_alpha)) / pow(mag_alpha, 2) * Alpha_ib_b.pow(2.);
	}
	else
	{
		C_new_old = Matrix3d::Identity() + Alpha_ib_b;
	}

	//SPECIFIC FORCE FRAME TRANSFORMATION
	Matrix3d ave_C_b_e;
	if (mag_alpha > 1e-8)
	{
		ave_C_b_e = m_INS_States.Cb_e*(Matrix3d::Identity() + (1 - cos(mag_alpha)) / pow(mag_alpha, 2)*Alpha_ib_b + (1 - sin(mag_alpha) / mag_alpha) / pow(mag_alpha, 2) * Alpha_ib_b.pow(2)) - 0.5*m_Skew(alpha_ie*Vector3d::UnitZ())*m_INS_States.Cb_e;
	}
	else
	{
		ave_C_b_e = m_INS_States.Cb_e - 0.5*m_Skew(alpha_ie*Vector3d::UnitZ())*m_INS_States.Cb_e;
	}

	Vector3d f_ib_e = ave_C_b_e*m_InsUserOutput->sf_actual;
	Vector3d gamma = GravityECEF(m_INS_States.position);

	// ATTITUDE UPDATE
	m_INS_States.Cb_e = C_Earth*m_INS_States.Cb_e*C_new_old;

	//POSITION UPDATE
	m_INS_States.position = m_INS_States.position + m_INS_States.velocity*dt + 0.5 * pow(dt, 2) * (f_ib_e + gamma - 2 * m_Skew(RtEarthRotn_Const*Vector3d::UnitZ()) * m_INS_States.velocity);

	// VELOCITY UPDATE
	m_INS_States.velocity = m_INS_States.velocity + (f_ib_e + gamma - 2 * m_Skew(RtEarthRotn_Const*Vector3d::UnitZ()) * m_INS_States.velocity)*dt;

	this->Calc_NED_States(m_INS_States);
}
/********************************************************************************************************
* Function               : Calc_NED_States
* Description            : Calculates NED States(lat,long,alt,vel_ned,rot_ned) from ECEF States(pos,vel,rot)
* Function Parameter     : None
* Return value           : None
* Specific library calls : None
* Functions called       : None
* Assumptions            : None
* Reference              : None
* Version History        :
* <1.1><Amanpreetsingh><02/03/2019>
***********************************************************************************************************/
void CInt::Calc_NED_States(CMain::INS_States& m_INS_States)
{
	Vector3d r = m_INS_States.position;
	// CALCULATING LONGITUDE
	m_INS_States.longitude = atan2(r[1], r[0]);

	double k1 = sqrt(1.0 - pow(Eccentricity_Const, 2.0))*abs(r[2]);
	double k2 = pow(Eccentricity_Const, 2.0)*EqtRadiusOfEarth_Const;
	double beta = Vector2d(r[0], r[1]).norm();
	double E = (k1 - k2) / beta;
	double F = (k1 + k2) / beta;
	double P = 4. / 3.*(E*F + 1.);
	double Q = 2.*(E*E - F*F);
	double D = pow(P, 3) + pow(Q, 2);
	double V = pow(sqrt(D) - Q, 1. / 3.) - pow(sqrt(D) + Q, 1. / 3.);
	double G = 0.5*(sqrt(E*E + V) + E);
	double T = sqrt(G*G + (F - V * G) / (2 * G - E)) - G;

	// CALCULATING LATITUDE
	double signlat = (r(2) > 0) - (r(2) < 0);
	m_INS_States.latitude = signlat* atan((1. - T*T) / (2 * T * sqrt(1 - pow(Eccentricity_Const, 2.0))));

	// CALCULATING HEIGHT
	double L_b = m_INS_States.latitude;
	m_INS_States.height = (beta - EqtRadiusOfEarth_Const * T) * cos(L_b) + (r(2) - signlat * EqtRadiusOfEarth_Const * sqrt(1 - pow(Eccentricity_Const, 2.0))) * sin(L_b);

	// ECEF to NED Transform Matrix
	Matrix3d Ce_n = m_TransformMatrix(L_b, m_INS_States.longitude);

	// TRANSFORMING VELOCITY TO NED
	m_INS_States.velocity_n = Ce_n*m_INS_States.velocity;

	// TRANSFORMING ATTITUDE TO NED
	m_INS_States.Cb_n = Ce_n*m_INS_States.Cb_e;
}

void CInt::InitAttitude(CMain::InsOutput *m_InsOutput, CMain::INS_States& m_INS_States){
	Matrix3d delta_Cb_n, est_Cb_n;
	double deg_to_rad = EIGEN_PI / 180.;
	delta_Cb_n = m_TransformMatrix(-0.05*deg_to_rad, 0.04*deg_to_rad, 1 * deg_to_rad);
	est_Cb_n = delta_Cb_n * m_InsOutput->cmat_bn;
	m_INS_States.Cb_e = m_TransformMatrix(m_INS_States.latitude, m_INS_States.longitude).transpose()*est_Cb_n;
}