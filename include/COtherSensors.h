#pragma once
#include "CSensor.h"
#include "GeographicLib/MagneticModel.hpp"
#include "GeographicLib/Geoid.hpp"
#include "stdatmos.h"

class COtherSensors : protected CSensor
{
public:
	COtherSensors();
	COtherSensors(string);
	~COtherSensors();
	struct MagReading
	{
		Eigen::Vector3d mag_strength;
		double mag_heading;
	};
	void m_Measure(CMain::TimeVar time, CMain::UserMot state);
	void m_ReadConfigFromFile(string config_file_name);
	MagReading m_RecentMagMeasurement;
	MagReading m_RecentMagTruth;
	double m_RecentPressMeasurement;
	double m_RecentPressTruth;
private:
	Eigen::Matrix3d C_mb;
	Eigen::Vector3d b_local_psd;
	Eigen::Vector3d b_hard;
	Eigen::Matrix3d M_soft;

	double p_noise_psd;
	const double r_eConst = 6.356766e6;

	inline double getNextCol4(ifstream& file)
	{
		string data;
		getline(file, data);
		stringstream ss(data);
		getline(ss, data, ',');
		getline(ss, data, ',');
		getline(ss, data, ',');
		getline(ss, data, ',');
		return stod(data);
	}

	GeographicLib::MagneticModel mag;
	GeographicLib::Geoid geo;
};