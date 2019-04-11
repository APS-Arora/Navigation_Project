#pragma once
#include "CSensor.h"
#include "c_main.h"
//#include "c_Integration.h"
#define _USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <random>
#include <chrono>
#include <time.h>
#include <cmath>
#include <fstream>
#include <vector>
#include <Eigen\Dense>

using namespace std;
using namespace Eigen;

class INS :
	public CSensor
{
	double m_step;
	bool m_turn;
	Vector3d sf_ib, wb_ib;
	Matrix3d Cb_e;
	double earth_rotn_rate = 0.00007292115, eccen = 0.0818191908425;
	//Updated INS Solution
public:
	INS();
	~INS();
	double NoiseGen(string, double);
	Vector3d GravityECEF(Vector3d position);
	void IMU_meas();
	void SensorRead();
	void main_function(CMain::UserMot m_UserMotion, CMain::InsOutput *);
	void INS_Estimate();
	void Calc_NED_States();
	struct CMain::UserMot m_InsUserMotion;
	struct CMain::InsOutput *m_InsUserOutput = new CMain::InsOutput;
	struct INS_States
	{
		Vector3d position, velocity, velocity_n, accel_bias,gyro_bias;
		Matrix3d Cb_e, Cb_n;
		double latitude, longitude, height;
	} m_INS_States;
	struct SensorData {
		SensorData() : acc_orientation(3), gyro_orientation(3), acc_misal(3), gyro_misal(3) {}
		vector<double> acc_orientation, gyro_orientation, acc_misal, gyro_misal;
		Matrix<double, 3, 3> acc_mis_mat, gyro_mis_mat;
		Vector3d acc_spd, acc_inrun_spd, acc_tc, acc_beta, gyro_spd, gyro_inrun_spd, gyro_tc, gyro_beta;
		Vector3d acc_fixbias, acc_fxbias, acc_fybias, acc_fzbias, acc_fxsqbias, acc_fysqbias, acc_fzsqbias, acc_fxfybias, acc_fyfzbias, acc_fzfxbias;
		Vector3d gyro_fixbias, gyro_fxbias, gyro_fybias, gyro_fzbias, gyro_fxsqbias, gyro_fysqbias, gyro_fzsqbias, gyro_fxfybias, gyro_fyfzbias, gyro_fzfxbias;
		Vector3d acc_fxscale, acc_fxsqscale, acc_fyscale, acc_fysqscale, acc_fzscale, acc_fzsqscale;
		Vector3d gyro_fixscale, gyro_fxscale, gyro_fxsqscale, gyro_fyscale, gyro_fysqscale, gyro_fzscale, gyro_fzsqscale;
	};
	SensorData m_Sensor;
};
