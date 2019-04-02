#pragma once
#include "c_main.h"

class CSensor{
public:
	virtual void m_Measure(CMain::TimeVar time, CMain::UserMot state) = 0;
	virtual void m_ReadConfigFromFile(string config_file) = 0;
	virtual void m_ReadConfigFromFile(ifstream config_file) = 0;

protected:
	Eigen::Matrix3d m_TransformMatrix(double lat, double longi);
	Eigen::Matrix3d m_TransformMatrix(double roll, double pitch, double yaw);
	Eigen::Matrix3d m_Skew(Eigen::Vector3d);
	double m_randn(void);

	default_random_engine m_rand;
	normal_distribution<double> m_gaussDist;
};