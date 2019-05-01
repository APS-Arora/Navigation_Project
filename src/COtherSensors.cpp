#include "COtherSensors.h"

COtherSensors::COtherSensors() : mag("wmm2015"), geo("egm2008-5")
{
	geo.CacheAll();
}


COtherSensors::~COtherSensors()
{
}

COtherSensors::COtherSensors(string config_file_name) : mag("wmm2015"), geo("egm2008-5")
{
	geo.CacheAll();
	this->m_ReadConfigFromFile(config_file_name);
}

void COtherSensors::m_ReadConfigFromFile(string config_file_name)
{
	ifstream config_file(config_file_name);
	string data;
	double psi = this->getNextCol4(config_file);
	double theta = this->getNextCol4(config_file);
	double phi = this->getNextCol4(config_file);
	this->C_mb = m_TransformMatrix(phi,theta,psi);

	this->b_hard[0] = this->getNextCol4(config_file);
	this->M_soft(0, 0) = this->getNextCol4(config_file);
	this->M_soft(1, 0) = this->getNextCol4(config_file);
	this->M_soft(2, 0) = this->getNextCol4(config_file);

	this->b_hard[1] = this->getNextCol4(config_file);
	this->M_soft(0, 1) = this->getNextCol4(config_file);
	this->M_soft(1, 1) = this->getNextCol4(config_file);
	this->M_soft(2, 1) = this->getNextCol4(config_file);

	this->b_hard[2] = this->getNextCol4(config_file);
	this->M_soft(0, 2) = this->getNextCol4(config_file);
	this->M_soft(1, 2) = this->getNextCol4(config_file);
	this->M_soft(2, 2) = this->getNextCol4(config_file);

	this->b_local_psd << this->getNextCol4(config_file),
						 this->getNextCol4(config_file),
						 this->getNextCol4(config_file);

	config_file.ignore(numeric_limits<streamsize>::max(), '\n');
	config_file.ignore(numeric_limits<streamsize>::max(), '\n');
	this->p_noise_psd = this->getNextCol4(config_file);
}

void COtherSensors::m_Measure(CMain::TimeVar time, CMain::UserMot state)
{
	Eigen::Matrix3d C_mn = C_mb*m_TransformMatrix(state.roll,state.pitch,state.yaw);
	double decyear = time.year + (time.days_in_year - 1. + time.time_of_day / 86400.) / (365. + ((time.year % 4) ? 0. : 1.));
	Eigen::Vector3d b_earth,b_local,b_heading;
	double Bx, By, Bz;
	mag(decyear, state.lat * 180. / EIGEN_PI, state.longi * 180. / EIGEN_PI, state.height, Bx, By, Bz);
	// ENU to NED Components Conversion
	b_earth << By, Bx, -Bz;
	b_earth /= 100.;
	b_local << m_gaussDist(m_rand), m_gaussDist(m_rand), m_gaussDist(m_rand);
	b_local = (b_local_psd.cwiseSqrt()).cwiseProduct(b_local);

	m_RecentMagMeasurement.mag_strength = b_hard + (Eigen::Matrix3d::Identity() + M_soft)*C_mn*(b_earth + b_local);
	m_RecentMagTruth.mag_strength = C_mn*b_earth;

	double yaw_nm = atan2(C_mn(0, 1), C_mn(0, 0))* 180. / EIGEN_PI;
	double declination, dummy;
	GeographicLib::MagneticModel::FieldComponents(Bx, By, Bz, dummy, dummy, declination, dummy);
	m_RecentMagTruth.mag_heading = yaw_nm - declination;
	b_heading = (C_mn*AngleAxisd(yaw_nm, Vector3d::UnitZ())).transpose()*m_RecentMagMeasurement.mag_strength;
	m_RecentMagMeasurement.mag_heading = atan2(-b_heading(1), b_heading(0)) * 180. / EIGEN_PI;

	double h_ortho = geo.ConvertHeight(state.lat * 180. / EIGEN_PI, state.longi * 180. / EIGEN_PI, state.height, GeographicLib::Geoid::ELLIPSOIDTOGEOID);
	double delta;
	SimpleAtmosphere(h_ortho/1000., dummy, delta, dummy);
	m_RecentPressTruth = delta*PZERO;
	m_RecentPressMeasurement = delta*PZERO + sqrt(p_noise_psd)*m_gaussDist(m_rand);
}