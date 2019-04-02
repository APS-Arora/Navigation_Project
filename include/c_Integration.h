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
	struct InsPosVel
	{
		Vector3d position, velocity;
		Matrix3d Cb_e;
	};
	InsPosVel *m_InsPosVel;
	CInt();
	~CInt();
	//LeastSquare();
	//GPSPosVel();
	//KalmanFilter();
	//MainFunction();
};

