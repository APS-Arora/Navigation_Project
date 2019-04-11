#pragma once
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
#include "ConstVar.h"
#include "INS.h"

using namespace std;
using namespace Eigen;

class CInt
{
	//GNSS Measurements: Sat Pos Vel, Range, Pseudorange and Rates
	//GNSS Config & Biases
	//Integrated Solution
	//State Prop Vector
	//Covariance Matrix
public:
	struct ErrorState
	{
		Vector3d position,
				 velocity,
				 attitude,
				 accel_bias,
				 gyro_bias;
		double   clock_offset,
				 clock_drift;
	} m_ErrorState;
	Matrix<double, 17, 17> ErrorCov;

	struct NoiseConfig
	{
		double  accel_PSD,
				gyro_PSD,
				accel_bias_PSD,
				gyro_bias_PSD,
				clock_offset_PSD,
				clock_drift_PSD;
	} m_NoiseConfig;

	CInt();
	~CInt();
	void m_predict(INS::INS_States, CMain::InsOutput);
	Vector3d GravityECEF(Vector3d position);
	void LsPosVel(CMain::SatData *, INS::INS_States&);
	//LeastSquare();
	//GPSPosVel();
	//KalmanFilter();
	//MainFunction();
protected:
	double dt;
	Matrix3d m_Skew(Vector3d);
	struct GpsYumaData
	{
		double  sat_id,
		health,
		eccent,
		time,
		a_o_i,
		rt_r_asc,
		sqrt_semi_maj,
		r_asc,
		arg_of_pge,
		m_ini,
		af0s,
		af1s,
		week,
		semi_maj;
	};
	GpsYumaData *mp_Yuma;
	int YumaRead();
};

