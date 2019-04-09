#include "CSensor.h"

Eigen::Matrix3d CSensor::m_TransformMatrix(double lat, double longi)
{
	// Use to Transform ECEF Vector to NED Vector (C^n_e)
	return(m_TransformMatrix(0, -EIGEN_PI / 2 - lat, longi));
}

Eigen::Matrix3d CSensor::m_TransformMatrix(double roll, double pitch, double yaw)
{
	return((AngleAxisd(-roll, Vector3d::UnitX())
		* AngleAxisd(-pitch, Vector3d::UnitY())
		* AngleAxisd(-yaw, Vector3d::UnitZ())).toRotationMatrix());
}

Matrix3d CSensor::m_Skew(Vector3d a){
	Matrix3d A;
	A << 0, -a(2), a(1),
		a(2), 0, -a(0),
		-a(1), a(0), 0;
	return A;
}

double CSensor::m_randn(void)
{
	return(m_gaussDist(m_rand));
}