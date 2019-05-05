/********************************************************************************************************
* Description            : Base Class implementing all functions required by all sensors
* Specific library calls : GeographicLib for Geoid Data and Magnetic Model
* Classes                : CSensor
* Assumptions            : None
* Reference              : 
* Version History        :
* <1.1><Amanpreetsingh><02/03/2019>
***********************************************************************************************************/
#include "CSensor.h"
/********************************************************************************************************
* Function               : m_TransformMatrix
* Description            : Transformation Matrix from ECEF to NED
* Function Parameter     : None
* Return value           : None
* Specific library calls : None
* Functions called       : None
* Assumptions            : None
* Reference              : None
* Version History        :
* <1.1><Amanpreetsingh><02/03/2019>
***********************************************************************************************************/
Eigen::Matrix3d CSensor::m_TransformMatrix(double lat, double longi)
{
	// Use to Transform ECEF Vector to NED Vector (C^n_e)
	return(m_TransformMatrix(0, -EIGEN_PI / 2 - lat, longi));
}
/********************************************************************************************************
* Function               : m_TransformMatrix
* Description            : Euler Angles to Rotation Matrix
* Function Parameter     : None
* Return value           : None
* Specific library calls : None
* Functions called       : None
* Assumptions            : None
* Reference              : None
* Version History        :
* <1.1><Amanpreetsingh><02/03/2019>
***********************************************************************************************************/
Eigen::Matrix3d CSensor::m_TransformMatrix(double roll, double pitch, double yaw)
{
	return((AngleAxisd(-roll, Vector3d::UnitX())
		* AngleAxisd(-pitch, Vector3d::UnitY())
		* AngleAxisd(-yaw, Vector3d::UnitZ())).toRotationMatrix());
}
/********************************************************************************************************
* Function               : m_Skew
* Description            : Calculates Cross Product Matrix from Vector
* Function Parameter     : None
* Return value           : None
* Specific library calls : None
* Functions called       : None
* Assumptions            : None
* Reference              : None
* Version History        :
* <1.1><Amanpreetsingh><02/03/2019>
***********************************************************************************************************/
Matrix3d CSensor::m_Skew(Vector3d a){
	Matrix3d A;
	A << 0, -a(2), a(1),
		a(2), 0, -a(0),
		-a(1), a(0), 0;
	return A;
}
/********************************************************************************************************
* Function               : m_randn
* Description            : Generates Samples from Gaussian Distribution with mean 0 and variance 1
* Function Parameter     : None
* Return value           : None
* Specific library calls : None
* Functions called       : None
* Assumptions            : None
* Reference              : None
* Version History        :
* <1.1><Amanpreetsingh><02/03/2019>
***********************************************************************************************************/
double CSensor::m_randn(void)
{
	return(m_gaussDist(m_rand));
}