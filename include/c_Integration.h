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
#include "CDelayCalc.h"

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
				clock_drift_PSD,
				pseudo_range_PSD,
				pseudo_range_rate_PSD;
	} m_NoiseConfig;

	CInt();
	~CInt();
	void m_predict(INS::INS_States, CMain::InsOutput);
	void m_correct(INS::INS_States&, CMain::GNSS_Measurement*);
	Vector3d GravityECEF(Vector3d position);
	void LsPosVel(CMain::SatData *, INS::INS_States&);
	void InitErrorCov();
	//LeastSquare();
	//GPSPosVel();
	//KalmanFilter();
	//MainFunction();
protected:
	double dt;
	Matrix3d m_Skew(Vector3d);
	int no_sat;
	CDelayCalc m_Delay;
	CMain::DelayCalcParam m_DelayParams;
};

