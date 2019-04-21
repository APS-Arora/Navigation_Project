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

class CInt : protected CSensor
{
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
	void run(CMain::INS_States&, CMain::GNSS_Measurement*, CMain::InsOutput);
	void m_predict(CMain::INS_States, CMain::InsOutput);
	void m_correct(CMain::INS_States&, CMain::GNSS_Measurement*);
	void LsPosVel(CMain::GNSS_Measurement*, CMain::INS_States&);
	void InitErrorCov();
	void INS_Estimate(CMain::INS_States&, CMain::InsOutput*);
	void Calc_NED_States(CMain::INS_States&);
	void InitAttitude(CMain::InsOutput*, CMain::INS_States&);

protected:
	double dt;
	Vector3d GravityECEF(Vector3d position);
	Matrix3d m_Skew(Vector3d);
	int no_sat;
	// Dependency on CDelayCalc to be removed in future
	CDelayCalc m_Delay;
	CMain::DelayCalcParam m_DelayParams;
	bool once;
};

