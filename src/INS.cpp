#define _USE_MATH_DEFINES
#include "INS.h"
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

using namespace std;
using namespace Eigen;

INS::INS()
{
	m_step = 0.1;
	m_turn = true;
	//cout << "Yo";
}


INS::~INS()
{
}



double INS::NoiseGen(string flag, double rms){

	double gravity = 9.81;
	if (flag == "gyrobias"){
		// M_PI represents pi (22/7) in _USE_MATH_DEFINES
		return(rms*M_PI / (180 * 3600)*randn());
	}

	if (flag == "accbias"){
		return(rms*(gravity / 1000)*randn());
	}

	if (flag == "scale"){
		return(rms / 1000000 * randn());
	}

	if (flag == "misal"){
		return(rms / 1000 * randn());
	}
}

Vector3d INS::GravityECEF(double pos[3]){
	Vector3d g, vec, gamma, position;
	position << pos[0], pos[1], pos[2];
	double R_0 = 6378137, mu = pow(3.986004418, 14), J_2 = 0.001082627;
	double z_scale, mag_r;
	mag_r = sqrt(position.transpose()*position);
	if (!mag_r){
		g << 0, 0, 0;
	}
	else {
		z_scale = 5 * pow(position[2] / mag_r, 2);
		vec << (1 - z_scale)*position[0], (1 - z_scale)*position[1], (3 - z_scale)*position[2];
		gamma = -mu / pow(mag_r, 3)*(position + 1.5*J_2*pow(R_0 / mag_r, 2)*vec);
		g[0] = gamma[0] + earth_rotn_rate*position[0];
		g[1] = gamma[1] + earth_rotn_rate*position[1];
		g[2] = gamma[2];
	}
	return g;
}


void INS::IMU_meas(void){

	Vector3d gravity, we_ie;
	Matrix3d cmat_en, cmat_bn, omega_e_ie;

	// Defining the constants to be used in the program
	we_ie << 0, 0, earth_rotn_rate;
	omega_e_ie << 0, -we_ie[2], we_ie[1],
		we_ie[2], 0, -we_ie[0],
		-we_ie[1], we_ie[0], 0;
	// Calculating the acceleration due to gravity using the latitude
	gravity << 0, 0, 9.7803253359*(1 + 0.001931853*pow(sin(m_InsUserMotion.lat), 2)) / pow(1 - pow(eccen*(sin(m_InsUserMotion.lat)), 2), 0.5);
	vector<double> ang_coord = { m_InsUserMotion.lat, m_InsUserMotion.longi };
	vector<double> orientation = { m_InsUserMotion.roll, m_InsUserMotion.pitch, m_InsUserMotion.yaw };
	// Calculating the ECEF to NED transformation matrix
	cmat_en = Transform_Matrix(ang_coord, "en");
	// Specific force in ECEF frame
	Vector3d sforce = Matrix<double, 3, 1>(m_InsUserMotion.user_acc) - cmat_en*gravity + 2 * omega_e_ie*Matrix<double, 3, 1>(m_InsUserMotion.user_vel);
	// Transforming specific force to NED to body frame
	cmat_bn = Transform_Matrix(orientation, "ab");
	sf_ib = cmat_bn*cmat_en.transpose()*sforce;
	wb_ib = Matrix<double, 3, 1>(m_InsUserMotion.user_ang) + cmat_bn*cmat_en.transpose()*we_ie;
}



void INS::SensorRead(){
	// ifstream object reading the sensor config file
	ifstream sensor("Sensor_Data.csv");
	string headers, waste; char comma;
	// Ensuring the files are open
	if (!sensor.is_open()){
		cout << "Error in reading configuration file" << endl;
	}
	//Reading the sensor file
	while (sensor.good()){
		getline(sensor, waste, ',');
		sensor >> m_Sensor.acc_orientation[0] >> comma >> m_Sensor.acc_orientation[1] >> comma >> m_Sensor.acc_orientation[2] >> comma >> m_Sensor.gyro_orientation[0] >> comma >> m_Sensor.gyro_orientation[1] >> comma >> m_Sensor.gyro_orientation[2] >> comma;
		sensor >> m_Sensor.acc_misal[0] >> comma >> m_Sensor.acc_misal[1] >> comma >> m_Sensor.acc_misal[2] >> comma >> m_Sensor.gyro_misal[0] >> comma >> m_Sensor.gyro_misal[1] >> comma >> m_Sensor.gyro_misal[2] >> comma;
		for (int i = 0; i < 3; i++){
			sensor >> m_Sensor.acc_inrun_spd[i] >> comma >> m_Sensor.acc_tc[i] >> comma >> m_Sensor.acc_fixbias[i] >> comma >> m_Sensor.acc_fxbias[i] >> comma >> m_Sensor.acc_fxsqbias[i] >> comma >> m_Sensor.acc_fxfybias[i] >> comma >> m_Sensor.acc_fybias[i] >> comma >> m_Sensor.acc_fysqbias[i] >> comma >> m_Sensor.acc_fyfzbias[i] >> comma >> m_Sensor.acc_fzbias[i] >> comma >> m_Sensor.acc_fzsqbias[i] >> comma >> m_Sensor.acc_fzfxbias[i] >> comma;
			sensor >> m_Sensor.acc_fxscale[i] >> comma >> m_Sensor.acc_fxsqscale[i] >> comma >> m_Sensor.acc_fyscale[i] >> comma >> m_Sensor.acc_fysqscale[i] >> comma >> m_Sensor.acc_fzscale[i] >> comma >> m_Sensor.acc_fzsqscale[i] >> comma >> m_Sensor.acc_spd[i] >> comma;
			sensor >> m_Sensor.gyro_inrun_spd[i] >> comma >> m_Sensor.gyro_tc[i] >> comma >> m_Sensor.gyro_fixbias[i] >> comma >> m_Sensor.gyro_fxbias[i] >> comma >> m_Sensor.gyro_fxsqbias[i] >> comma >> m_Sensor.gyro_fxfybias[i] >> comma >> m_Sensor.gyro_fybias[i] >> comma >> m_Sensor.gyro_fysqbias[i] >> comma >> m_Sensor.gyro_fyfzbias[i] >> comma >> m_Sensor.gyro_fzbias[i] >> comma >> m_Sensor.gyro_fzsqbias[i] >> comma >> m_Sensor.gyro_fzfxbias[i] >> comma;
			sensor >> m_Sensor.gyro_fixscale[i] >> comma >> m_Sensor.gyro_fxscale[i] >> comma >> m_Sensor.gyro_fxsqscale[i] >> comma >> m_Sensor.gyro_fyscale[i] >> comma >> m_Sensor.gyro_fysqscale[i] >> comma >> m_Sensor.gyro_fzscale[i] >> comma >> m_Sensor.gyro_fzsqscale[i] >> comma >> m_Sensor.gyro_spd[i] >> comma;
		}
	}
	sensor.close();

	//Computing beta and spectral power density for Gauss-Markov process to be used for computing in-run bias
	m_Sensor.acc_spd = m_Sensor.acc_spd / (9.81 * 1000);
	m_Sensor.gyro_spd = m_Sensor.gyro_spd / (60 * 180 / M_PI);
	for (int j = 0; j < 3; j++){
		m_Sensor.acc_beta[j] = 1 / m_Sensor.acc_tc[j];
		m_Sensor.gyro_beta[j] = 1 / m_Sensor.gyro_tc[j];
		m_Sensor.acc_inrun_spd[j] = (m_Sensor.acc_inrun_spd[j] * (1 - m_Sensor.acc_beta[j] * m_step / 2) / m_step) / (9.81 * 1000);
		m_Sensor.gyro_inrun_spd[j] = (m_Sensor.gyro_inrun_spd[j] * (1 - m_Sensor.gyro_beta[j] * m_step / 2) / m_step) / (60 * 180 / M_PI);
	}
	// ACC : Inititalizing the error values from commercial sensor data using NoiseGen
	vector<double> acc_am = { NoiseGen("misal", m_Sensor.acc_misal[0]), NoiseGen("misal", m_Sensor.acc_misal[1]), NoiseGen("misal", m_Sensor.acc_misal[2]) };
	m_Sensor.acc_fixbias = { NoiseGen("accbias", m_Sensor.acc_fixbias[0]), NoiseGen("accbias", m_Sensor.acc_fixbias[1]), NoiseGen("accbias", m_Sensor.acc_fixbias[2]) };
	m_Sensor.acc_fxbias = { NoiseGen("accbias", m_Sensor.acc_fxbias[0]), NoiseGen("accbias", m_Sensor.acc_fxbias[1]), NoiseGen("accbias", m_Sensor.acc_fxbias[2]) };
	m_Sensor.acc_fybias = { NoiseGen("accbias", m_Sensor.acc_fybias[0]), NoiseGen("accbias", m_Sensor.acc_fybias[1]), NoiseGen("accbias", m_Sensor.acc_fybias[2]) };
	m_Sensor.acc_fzbias = { NoiseGen("accbias", m_Sensor.acc_fzbias[0]), NoiseGen("accbias", m_Sensor.acc_fzbias[1]), NoiseGen("accbias", m_Sensor.acc_fzbias[2]) };
	m_Sensor.acc_fxsqbias = { NoiseGen("accbias", m_Sensor.acc_fxsqbias[0]), NoiseGen("accbias", m_Sensor.acc_fxsqbias[1]), NoiseGen("accbias", m_Sensor.acc_fxsqbias[2]) };
	m_Sensor.acc_fysqbias = { NoiseGen("accbias", m_Sensor.acc_fysqbias[0]), NoiseGen("accbias", m_Sensor.acc_fysqbias[1]), NoiseGen("accbias", m_Sensor.acc_fysqbias[2]) };
	m_Sensor.acc_fzsqbias = { NoiseGen("accbias", m_Sensor.acc_fzsqbias[0]), NoiseGen("accbias", m_Sensor.acc_fzsqbias[1]), NoiseGen("accbias", m_Sensor.acc_fzsqbias[2]) };
	m_Sensor.acc_fxfybias = { NoiseGen("accbias", m_Sensor.acc_fxfybias[0]), NoiseGen("accbias", m_Sensor.acc_fxfybias[1]), NoiseGen("accbias", m_Sensor.acc_fxfybias[2]) };
	m_Sensor.acc_fyfzbias = { NoiseGen("accbias", m_Sensor.acc_fyfzbias[0]), NoiseGen("accbias", m_Sensor.acc_fyfzbias[1]), NoiseGen("accbias", m_Sensor.acc_fyfzbias[2]) };
	m_Sensor.acc_fzfxbias = { NoiseGen("accbias", m_Sensor.acc_fzfxbias[0]), NoiseGen("accbias", m_Sensor.acc_fzfxbias[1]), NoiseGen("accbias", m_Sensor.acc_fzfxbias[2]) };
	m_Sensor.acc_fxscale = { NoiseGen("scale", m_Sensor.acc_fxscale[0]), NoiseGen("scale", m_Sensor.acc_fxscale[1]), NoiseGen("scale", m_Sensor.acc_fxscale[2]) };
	m_Sensor.acc_fyscale = { NoiseGen("scale", m_Sensor.acc_fyscale[0]), NoiseGen("scale", m_Sensor.acc_fyscale[1]), NoiseGen("scale", m_Sensor.acc_fyscale[2]) };
	m_Sensor.acc_fzscale = { NoiseGen("scale", m_Sensor.acc_fzscale[0]), NoiseGen("scale", m_Sensor.acc_fzscale[1]), NoiseGen("scale", m_Sensor.acc_fzscale[2]) };
	m_Sensor.acc_fxsqscale = { NoiseGen("scale", m_Sensor.acc_fxsqscale[0]), NoiseGen("scale", m_Sensor.acc_fxsqscale[1]), NoiseGen("scale", m_Sensor.acc_fxsqscale[2]) };
	m_Sensor.acc_fysqscale = { NoiseGen("scale", m_Sensor.acc_fysqscale[0]), NoiseGen("scale", m_Sensor.acc_fysqscale[1]), NoiseGen("scale", m_Sensor.acc_fysqscale[2]) };
	m_Sensor.acc_fzsqscale = { NoiseGen("scale", m_Sensor.acc_fzsqscale[0]), NoiseGen("scale", m_Sensor.acc_fzsqscale[1]), NoiseGen("scale", m_Sensor.acc_fzsqscale[2]) };
	m_Sensor.acc_mis_mat = Transform_Matrix(acc_am, "ab");

	// GYRO : Inititalizing the error values from commercial sensor data using NoiseGen
	vector<double> gyro_am = { NoiseGen("misal", m_Sensor.gyro_misal[0]), NoiseGen("misal", m_Sensor.gyro_misal[1]), NoiseGen("misal", m_Sensor.gyro_misal[2]) };
	m_Sensor.gyro_fixbias = { NoiseGen("gyrobias", m_Sensor.gyro_fixbias[0]), NoiseGen("gyrobias", m_Sensor.gyro_fixbias[1]), NoiseGen("gyrobias", m_Sensor.gyro_fixbias[2]) };
	m_Sensor.gyro_fxbias = { NoiseGen("gyrobias", m_Sensor.gyro_fxbias[0]), NoiseGen("gyrobias", m_Sensor.gyro_fxbias[1]), NoiseGen("gyrobias", m_Sensor.gyro_fxbias[2]) };
	m_Sensor.gyro_fybias = { NoiseGen("gyrobias", m_Sensor.gyro_fybias[0]), NoiseGen("gyrobias", m_Sensor.gyro_fybias[1]), NoiseGen("gyrobias", m_Sensor.gyro_fybias[2]) };
	m_Sensor.gyro_fzbias = { NoiseGen("gyrobias", m_Sensor.gyro_fzbias[0]), NoiseGen("gyrobias", m_Sensor.gyro_fzbias[1]), NoiseGen("gyrobias", m_Sensor.gyro_fzbias[2]) };
	m_Sensor.gyro_fxsqbias = { NoiseGen("gyrobias", m_Sensor.gyro_fxsqbias[0]), NoiseGen("gyrobias", m_Sensor.gyro_fxsqbias[1]), NoiseGen("gyrobias", m_Sensor.gyro_fxsqbias[2]) };
	m_Sensor.gyro_fysqbias = { NoiseGen("gyrobias", m_Sensor.gyro_fysqbias[0]), NoiseGen("gyrobias", m_Sensor.gyro_fysqbias[1]), NoiseGen("gyrobias", m_Sensor.gyro_fysqbias[2]) };
	m_Sensor.gyro_fzsqbias = { NoiseGen("gyrobias", m_Sensor.gyro_fzsqbias[0]), NoiseGen("gyrobias", m_Sensor.gyro_fzsqbias[1]), NoiseGen("gyrobias", m_Sensor.gyro_fzsqbias[2]) };
	m_Sensor.gyro_fxfybias = { NoiseGen("gyrobias", m_Sensor.gyro_fxfybias[0]), NoiseGen("gyrobias", m_Sensor.gyro_fxfybias[1]), NoiseGen("gyrobias", m_Sensor.gyro_fxfybias[2]) };
	m_Sensor.gyro_fyfzbias = { NoiseGen("gyrobias", m_Sensor.gyro_fyfzbias[0]), NoiseGen("gyrobias", m_Sensor.gyro_fyfzbias[1]), NoiseGen("gyrobias", m_Sensor.gyro_fyfzbias[2]) };
	m_Sensor.gyro_fzfxbias = { NoiseGen("gyrobias", m_Sensor.gyro_fzfxbias[0]), NoiseGen("gyrobias", m_Sensor.gyro_fzfxbias[1]), NoiseGen("gyrobias", m_Sensor.gyro_fzfxbias[2]) };
	m_Sensor.gyro_fixscale = { NoiseGen("scale", m_Sensor.gyro_fixscale[0]), NoiseGen("scale", m_Sensor.gyro_fixscale[1]), NoiseGen("scale", m_Sensor.gyro_fixscale[2]) };
	m_Sensor.gyro_fxscale = { NoiseGen("scale", m_Sensor.gyro_fxscale[0]), NoiseGen("scale", m_Sensor.gyro_fxscale[1]), NoiseGen("scale", m_Sensor.gyro_fxscale[2]) };
	m_Sensor.gyro_fyscale = { NoiseGen("scale", m_Sensor.gyro_fyscale[0]), NoiseGen("scale", m_Sensor.gyro_fyscale[1]), NoiseGen("scale", m_Sensor.gyro_fyscale[2]) };
	m_Sensor.gyro_fzscale = { NoiseGen("scale", m_Sensor.gyro_fzscale[0]), NoiseGen("scale", m_Sensor.gyro_fzscale[1]), NoiseGen("scale", m_Sensor.gyro_fzscale[2]) };
	m_Sensor.gyro_fxsqscale = { NoiseGen("scale", m_Sensor.gyro_fxsqscale[0]), NoiseGen("scale", m_Sensor.gyro_fxsqscale[1]), NoiseGen("scale", m_Sensor.gyro_fxsqscale[2]) };
	m_Sensor.gyro_fysqscale = { NoiseGen("scale", m_Sensor.gyro_fysqscale[0]), NoiseGen("scale", m_Sensor.gyro_fysqscale[1]), NoiseGen("scale", m_Sensor.gyro_fysqscale[2]) };
	m_Sensor.gyro_fzsqscale = { NoiseGen("scale", m_Sensor.gyro_fzsqscale[0]), NoiseGen("scale", m_Sensor.gyro_fzsqscale[1]), NoiseGen("scale", m_Sensor.gyro_fzsqscale[2]) };
	m_Sensor.gyro_mis_mat = Transform_Matrix(gyro_am, "ab");

	// Aborting calling of this function after the first iteration
	m_InsUserOutput->acc_inrun = { 0, 0, 0 }, m_InsUserOutput->gyro_inrun = { 0, 0, 0 };
	m_turn = false;
	//cout << "Success";
}



void INS::main_function(CMain::UserMot m_UserMotion, CMain::InsOutput *m_InsOutput){

	memcpy(&m_InsUserMotion, &m_UserMotion, sizeof(m_InsUserMotion));

	// Calling function for calculating truth
	IMU_meas();

	//Function call to read sensor file
	if (m_turn == true)
	{
		SensorRead();
	}
	//Local variables for storing data
		Vector3d sf_acc = Transform_Matrix(m_Sensor.acc_orientation, "ab")*sf_ib;
		Vector3d ang_gyro = Transform_Matrix(m_Sensor.gyro_orientation, "ab")*wb_ib;
		Vector3d gyro_scale;

		//Computing the errors
		for (size_t j = 0; j < 3; j++){
			m_InsUserOutput->acc_bias[j] = m_Sensor.acc_fixbias[j] + m_Sensor.acc_fxbias[j] * sf_acc[0] + m_Sensor.acc_fybias[j] * sf_acc[1]
				+ m_Sensor.acc_fzbias[j] * sf_acc[2] + m_Sensor.acc_fxsqbias[j] * (sf_acc[0] * sf_acc[0]) + m_Sensor.acc_fysqbias[j] * (sf_acc[1] * sf_acc[1])
				+ m_Sensor.acc_fzsqbias[j] * (sf_acc[2] * sf_acc[2]) + m_Sensor.acc_fxfybias[j] * (sf_acc[0] * sf_acc[1]) +
				m_Sensor.acc_fyfzbias[j] * (sf_acc[1] * sf_acc[2]) + m_Sensor.acc_fzfxbias[j] * (sf_acc[2] * sf_acc[0]);
			m_InsUserOutput->acc_scale[j] = m_Sensor.acc_fxscale[j] * sf_acc[0] + m_Sensor.acc_fyscale[j] * sf_acc[1] + m_Sensor.acc_fzscale[j] * sf_acc[2]
				+ m_Sensor.acc_fxsqscale[j] * (sf_acc[0] * sf_acc[0]) + m_Sensor.acc_fysqscale[j] * (sf_acc[1] * sf_acc[1]) + m_Sensor.acc_fzsqscale[j] * (sf_acc[2] * sf_acc[2]);
			m_InsUserOutput->gyro_bias[j] = m_Sensor.gyro_fixbias[j] + m_Sensor.gyro_fxbias[j] * sf_acc[0] + m_Sensor.gyro_fybias[j] * sf_acc[1]
				+ m_Sensor.gyro_fzbias[j] * sf_acc[2] + m_Sensor.gyro_fxsqbias[j] * (sf_acc[0] * sf_acc[0]) + m_Sensor.gyro_fysqbias[j] * (sf_acc[1] * sf_acc[1])
				+ m_Sensor.gyro_fzsqbias[j] * (sf_acc[2] * sf_acc[2]) + m_Sensor.gyro_fxfybias[j] * (sf_acc[0] * sf_acc[1]) +
				m_Sensor.gyro_fyfzbias[j] * (sf_acc[1] * sf_acc[2]) + m_Sensor.gyro_fzfxbias[j] * (sf_acc[2] * sf_acc[0]);
			gyro_scale[j] = m_Sensor.gyro_fixscale[j] + m_Sensor.gyro_fxscale[j] * sf_acc[0] + m_Sensor.gyro_fyscale[j] * sf_acc[1] + m_Sensor.gyro_fzscale[j] * sf_acc[2]
				+ m_Sensor.gyro_fxsqscale[j] * (sf_acc[0] * sf_acc[0]) + m_Sensor.gyro_fysqscale[j] * (sf_acc[1] * sf_acc[1]) + m_Sensor.gyro_fzsqscale[j] * (sf_acc[2] * sf_acc[2]);
			m_InsUserOutput->acc_inrun[j] = (1 - (m_Sensor.acc_beta[j] * m_step))*m_InsUserOutput->acc_inrun[j] + sqrt(2 * m_Sensor.acc_beta[j])*(m_Sensor.acc_inrun_spd[j] * randn()*m_step);
			m_InsUserOutput->gyro_inrun[j] = (1 - (m_Sensor.gyro_beta[j] * m_step))*m_InsUserOutput->gyro_inrun[j] + sqrt(2 * m_Sensor.gyro_beta[j])*(m_Sensor.gyro_inrun_spd[j] * randn()*m_step);
			m_InsUserOutput->acc_noise[j] = m_Sensor.acc_spd[j] * randn()*m_step;
			m_InsUserOutput->gyro_noise[j] = m_Sensor.gyro_spd[j] * randn()*m_step;
		}
		
		//Isolating the axis misalignment
		for (size_t j = 0; j < 3; j++)
		{
			m_InsUserOutput->acc_misal_error[j] = (m_Sensor.acc_mis_mat*sf_acc - sf_acc)[j];
		}

		Matrix<double, 3, 3> gyro_scale_diag;
		gyro_scale_diag << gyro_scale[0], 0, 0, 0, gyro_scale[1], 0, 0, 0, gyro_scale[2];
		Matrix3d gyro_error_mat = m_Sensor.gyro_mis_mat + Matrix3d::Identity() + gyro_scale_diag;

		//Computing gyro scale factor and misaligment error
		m_InsUserOutput->gyro_scale_error = gyro_scale_diag*ang_gyro;
		m_InsUserOutput->gyro_misal_error = m_Sensor.gyro_mis_mat*ang_gyro - ang_gyro;
		// Adding the error terms to the truth
		m_InsUserOutput->sf_actual = sf_acc + (m_InsUserOutput->acc_bias + m_InsUserOutput->acc_misal_error + m_InsUserOutput->acc_inrun + m_InsUserOutput->acc_scale + m_InsUserOutput->acc_noise);
		m_InsUserOutput->ang_actual = ang_gyro + m_InsUserOutput->gyro_bias + gyro_error_mat*ang_gyro + m_InsUserOutput->gyro_inrun + m_InsUserOutput->gyro_noise;

		
		//Copying the data to Cmain object
		memcpy(m_InsOutput, m_InsUserOutput, sizeof(*m_InsOutput));
		//return(m_InsUserOutput);
}

void INS::EstPosVel(struct CMain::UserMot m_UserMotion, struct CInt::InsPosVel *m_InsPosVel){

	double alpha_ie = earth_rotn_rate * m_step, mag_alpha ;
	Vector3d alpha_ib_b, alpha_ie_vec, earth_vec, f_ib_e, user_vel, user_pos;
	user_vel << m_UserMotion.user_vel[0], m_UserMotion.user_vel[1], m_UserMotion.user_vel[2];
	user_pos << m_UserMotion.user_pos[0], m_UserMotion.user_pos[1], m_UserMotion.user_pos[2];
	alpha_ie_vec << 0, 0, alpha_ie;
	earth_vec << 0, 0, earth_rotn_rate;
	Matrix3d C_Earth, Alpha_ib_b, Cnew_old, ave_Cb_e;

	C_Earth << cos(alpha_ie), sin(alpha_ie), 0,
		- sin(alpha_ie), cos(alpha_ie), 0,
		0, 0, 1;
	alpha_ib_b = m_InsUserOutput->ang_actual * m_step;
	mag_alpha = sqrt(alpha_ib_b.transpose() * alpha_ib_b);
	Alpha_ib_b = Skew(alpha_ib_b);

	if (mag_alpha > 0.00000001){
		Cnew_old = Matrix3d::Identity() + sin(mag_alpha) / mag_alpha * Alpha_ib_b + (1 - cos(mag_alpha)) / mag_alpha *mag_alpha * Alpha_ib_b * Alpha_ib_b;
	}
	else Cnew_old = Matrix3d::Identity() + Alpha_ib_b;

	m_InsUserPosVel->Cb_e = C_Earth * Cb_e * Cnew_old;

	if (mag_alpha > 0.00000001){
		ave_Cb_e = Cb_e * (Matrix3d::Identity() + (1 - cos(mag_alpha)) / mag_alpha *mag_alpha
			* Alpha_ib_b + (1 - sin(mag_alpha) / mag_alpha) / mag_alpha *mag_alpha
			* Alpha_ib_b * Alpha_ib_b) - 0.5 * Skew(alpha_ie_vec)
			* Cb_e;
	}
	else ave_Cb_e = Cb_e - 0.5 * Skew(alpha_ie_vec) *Cb_e;
	//Transform specific force to ECEF - frame resolving axes using (5.85)
	f_ib_e = ave_Cb_e * m_InsUserOutput->sf_actual;
	//Updating Velocity
	m_InsUserPosVel->velocity = user_vel + m_step * (f_ib_e + GravityECEF(m_UserMotion.user_pos) - 2 * Skew(earth_vec) * user_vel);
	//Updating Position
	m_InsUserPosVel->position = user_pos + (m_InsUserPosVel->velocity + user_vel)*0.5*m_step;
}