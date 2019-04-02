#include "Sensor.h"
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

Sensor::Sensor()
{
}


Sensor::~Sensor()
{
}

// Function for producing random values using standard normal distribution
double Sensor::randn(void){
	//double value;
	//normal_distribution<double> x(0, 1);
	//unsigned seed = chrono::steady_clock::now().time_since_epoch().count();
	//default_random_engine eng(seed);
	//value = x(eng);
	return(2.04312);
}

// Function for calculating the transformation matrix with given vector of rotations
Matrix3d Sensor::Transform_Matrix(vector<double> V, string flag){
	Matrix<double, 3, 3> trans_mat;
	if (flag == "en"){
		trans_mat << -sin(V[0])*cos(V[1]), -sin(V[1]), cos(V[0])*cos(V[1]),
			-sin(V[0])*sin(V[1]), cos(V[1]), -cos(V[0])*sin(V[1]),
			cos(V[0]), 0, -sin(V[0]);
	}
	if (flag == "ab"){
		trans_mat << cos(V[1])*cos(V[2]), cos(V[1])*sin(V[2]), -sin(V[1]),
			-cos(V[0])*sin(V[2]) + sin(V[0])*sin(V[1])*cos(V[2]), cos(V[0])*cos(V[2]) + sin(V[0])*sin(V[1])*sin(V[2]), sin(V[0])*cos(V[1]),
			sin(V[0])*sin(V[2]) + cos(V[0])*sin(V[1])*cos(V[2]), -sin(V[0])*cos(V[2]) + cos(V[0])*sin(V[1])*sin(V[2]), cos(V[0])*cos(V[1]);
	}
	return(trans_mat);
}

Matrix3d Sensor::Skew(Vector3d a){
	Matrix3d A;
	A << 0, -a(3), a(2),
		a(3), 0, -a(1),
	- a(2), a(1), 0;
	return A;
}