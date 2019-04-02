#pragma once
#include "c_main.h"
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
// Header file for the Sensor class
using namespace Eigen;

#ifndef Sensor_H
#define Sensor_H
class Sensor
{
public:
	Sensor();
	~Sensor();
	double randn(void);
	Matrix3d Skew(Vector3d);
	Matrix3d Transform_Matrix(vector<double>, string);
};

#endif

