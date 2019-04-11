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

int CInt::YumaRead()
{

	// take input from file 
	ifstream iFile("gps_yuma.txt");

	// Output file 
	ofstream ofile_1("gps_orbital.txt");

	string strn;
	double inte;

	// Commands to read the Yuma file, reading word by word in each line 
	int count = 0;
	int index = -1;
	while (!iFile.eof())
	{
		string temp;
		iFile >> temp;
		if (iFile.eof()) break;

		if (temp == "********")
		{
			iFile >> strn; iFile >> strn; iFile >> strn; iFile >> strn; iFile >> strn; iFile >> strn;
			index++;
			count = 0;
			continue;
		}
		if (temp == "ID:")
		{
			iFile >> inte;
			mp_Yuma[index].sat_id = inte;
		}
		else if (temp == "Health:")
		{
			iFile >> inte;
			mp_Yuma[index].health = inte;
		}
		else if (temp == "Eccentricity:")
		{
			iFile >> inte;
			mp_Yuma[index].eccent = inte;
		}
		else if (temp == "Time")
		{
			iFile >> strn; iFile >> strn; iFile >> inte;
			mp_Yuma[index].time = inte;
		}
		else if (temp == "Orbital")
		{
			iFile >> strn; iFile >> inte;
			mp_Yuma[index].a_o_i = inte;
		}
		else if (temp == "Rate")
		{
			iFile >> strn; iFile >> strn; iFile >> strn; iFile >> inte;
			mp_Yuma[index].rt_r_asc = inte;
		}
		else if (temp == "SQRT(A)")
		{
			iFile >> strn; iFile >> strn; iFile >> inte;

			mp_Yuma[index].sqrt_semi_maj = inte;
			mp_Yuma[index].semi_maj = (inte*inte);
		}
		else if (temp == "Right")
		{
			iFile >> strn; iFile >> strn; iFile >> strn; iFile >> inte;
			mp_Yuma[index].r_asc = inte;
		}
		else if (temp == "Argument")
		{
			iFile >> strn; iFile >> strn; iFile >> inte;
			mp_Yuma[index].arg_of_pge = inte;
		}
		else if (temp == "Mean")
		{
			iFile >> strn;
			iFile >> inte;
			mp_Yuma[index].m_ini = inte;
		}
		else if (temp == "Af0(s):")
		{
			iFile >> inte;
			mp_Yuma[index].af0s = inte;
		}
		else if (temp == "Af1(s/s):")
		{

			iFile >> inte;
			mp_Yuma[index].af1s = inte;
		}
		else if (temp == "week:")
		{

			iFile >> inte;
			mp_Yuma[index].week = inte;
		}

	}
	// Save the yuma file data in a file 
	for (int iter = 0; iter <= index; iter++)
	{
		ofile_1 << mp_Yuma[iter].sat_id << "  " << mp_Yuma[iter].health << "  " << mp_Yuma[iter].eccent
			<< "  " << mp_Yuma[iter].time << "  " << mp_Yuma[iter].a_o_i << "  " << mp_Yuma[iter].rt_r_asc
			<< "  " << mp_Yuma[iter].semi_maj << "  " << mp_Yuma[iter].r_asc << "  " << mp_Yuma[iter].arg_of_pge
			<< "  " << mp_Yuma[iter].m_ini << "  " << mp_Yuma[iter].af0s << "  " << mp_Yuma[iter].af1s << "  "
			<< mp_Yuma[iter].week << endl;
		ofile_1 << "#$$$$$\n";

	}
	return index;
}
/*
void CInt::LsPosVel(CMain::SatData *mp_GpsSatData, INS::INS_States& INS_Init){
	//Calculating number of satellites
	int no_sat = YumaRead();
	//Speed of Light
	int c = 299792458;
	int omega_ie = 0.00007292115;
	Vector4d x_pred = { 0, 0, 0, 0 }, x_est;
	Matrix<double, Dynamic, 4> H_mat = Matrix<double, Dynamic, 4>::Zero(no_sat,4);
	Matrix<double, Dynamic, 3> gnss_pos = Matrix<double, Dynamic, 3>::Zero(no_sat, 3), 
							   gnss_vel = Matrix<double, Dynamic, 3>::Zero(no_sat, 3);
	VectorXd pred_meas = VectorXd::Zero(no_sat), gnss_range = VectorXd::Zero(no_sat), gnss_range_rate = VectorXd::Zero(no_sat);
	Vector3d delta_r, est_pos, omega, u_as_e;
	Matrix3d Ce_i, Omega_ie;
	double approx_range, range, range_rate;
	int cnvg = 1;

	while (cnvg > 0.0001){
		for (int j = 0; j < no_sat; j++){
			gnss_pos(j, 0) = mp_GpsSatData[j].x_cord; gnss_pos(j, 1) = mp_GpsSatData[j].y_cord, gnss_pos(j, 2) = mp_GpsSatData[j].z_cord;
			gnss_range(j) = mp_GpsSatData[j].range;
			//Predict approx range
			delta_r = gnss_pos.row(j) - x_pred.head(3);
			approx_range = sqrt(delta_r.transpose()*delta_r);

			//Frame rotation during signal transit time (8.36)
			Ce_i << 1, omega_ie * approx_range / c, 0,
				-omega_ie * approx_range / c, 1, 0,
				0, 0, 1;

			//Predict pseudorange (9.144)
			delta_r = Ce_i*gnss_pos.row(j) - x_pred.head(3);
			range = delta_r.norm();
			pred_meas[j] = range + x_pred(4);

			//Predict line of sight and deploy in measurement matrix, (9.144)
			H_mat.block<1, 3>(j, 0) = -delta_r.transpose() / range;
			H_mat(j, 4) = 1;
		}
		//Unweighted least - squares solution
		x_est = x_pred + (H_mat.transpose()*H_mat).inverse()*H_mat.transpose()*(gnss_range - pred_meas);

		//Test convergence
		cnvg = sqrt((x_est - x_pred).transpose()*(x_est - x_pred));

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
			gnss_vel(j, 0) = mp_GpsSatData[j].x_vel; gnss_vel(j, 0) = mp_GpsSatData[j].y_vel; gnss_vel(j, 0) = mp_GpsSatData[j].z_vel;
			gnss_range_rate(j) = mp_GpsSatData[j].range_rt;
			//Predict approx range
			delta_r = gnss_pos.row(j) - x_pred.head(3);
			approx_range = sqrt(delta_r.transpose()*delta_r);

			//Frame rotation during signal transit time (8.36)
			Ce_i << 1, omega_ie * approx_range / c, 0,
				-omega_ie * approx_range / c, 1, 0,
				0, 0, 1;

			//Predict pseudorange (9.144)
			delta_r = Ce_i*gnss_pos.row(j) - x_pred.head(3);
			range = delta_r.norm();

			//Calculate line of sight using (8.41)
			u_as_e = delta_r / range;

			//Predict range rate
			range_rate = u_as_e.transpose()*(Ce_i*gnss_vel.row(j).transpose() + omega_ie*gnss_pos.row(j).transpose() - x_pred.head(3) + Omega_ie*INS_Init.position);
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
*/