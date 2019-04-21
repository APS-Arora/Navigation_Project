/********************************************************************************************************
* File				     : CDelayCalc.cpp
*
* Description            : This cpp file computes Atmospheric Delays
*							1. Ionospheric Delay: Klobuchar Model
*							2. Tropospheric Delay: RTCA98 and Hopfield model
*
* Specific library calls : None
*
* Classes                : CDelayCalc
*
* Functions called       : ComputeDelays()
*						   TropoDelayCalcRTCA()
*						   TropoDelayRTCA()
*						   IonoDelayKlob()
*						   TropoDelayHop()
*						   ComputeWeight()
*						   AccessTable()
*						   ComputeMeanSeaCorrection()
*						   StdAtm()
*						 
* Assumptions            : 1. Height of antenna above sea level is at max 30km 
*						   2. Invalid Conditions are not specified in the code. It is assumed that
*						      user will give proper inputs
*
* Reference              : 1. Klobuchar Model: GPS Interface Control Document IS-GPS-200H
*						   2. Hopfield Model: GPS Theory and Practice (3rd revised edition) by 
*						      B.Hofmann-Wellenhof, H.Lichtenegger, and J.Collins 
*						   3. RTCA98 Model: 
*
* Version History        :
* <1.1><Thakur Shivam Singh><29/9/2017>
********************************************************************************************************/

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<math.h>
#include<stdlib.h>
#include<Eigen/Dense>
#include "c_Gps.h"
#include "c_Irnss.h"
#include "ConstVar.h"
#include "c_main.h"
#include "CDelayCalc.h"
#include <vector>


using namespace std;
using namespace Eigen;

CDelayCalc::CDelayCalc() : geo("egm2008-5")
{
}

CDelayCalc::~CDelayCalc()
{
}

/********************************************************************************************************
* Description            : Function to calculate Atmospheric Delays
*						   This is a public member function which is called from outside the class
*
* Specific library calls : None
*
* Classes                : CDelayCalc
*
* Functions called       : TropoDelayRTCA()
*                          IonoDelayKlob()
*						   TropoDelayHop()
*
* Assumptions            : None
*
* Version History        :
* <1.1><Thakur Shivam Singh><29/9/2017>
********************************************************************************************************/

void CDelayCalc::ComputeDelays(double time_of_week, 
								double day_of_year, 
								double lati, 
								double longi, 
								double height,
								double azimuth,
								double elevation, 
								struct AtmSatDelay &m_DelayCalc,
								struct DelayCalcParam &DelayParam)
{
	
	double mask_angle;
	// Copy the values obtained from the function parameters
	m_Lat = lati;
	m_Longi = longi;
	m_Ele = elevation;
	m_Azi = azimuth;
	m_Height = height;
	m_TimeOfWeek = time_of_week;
	m_DayOfYear = day_of_year;

	for (int iter = 0; iter < 4; iter++)
	{
		m_Alpha[iter] = DelayParam.alpha[iter];
		m_Beta[iter] = DelayParam.beta[iter];
	}
	m_Humidity = DelayParam.humidity;
	mask_angle = DelayParam.mask_angle * DegreetoRadian_Const;

	// Atmospheric delays are not calculated if elevation angle is less than mask angle
	if (m_Ele > mask_angle)
	{
		TropoDelayHop();
		IonoDelayKlob();
	}
	else
	{
		m_AtmDelayCalc.tropo_delay_hop = 0;
		m_AtmDelayCalc.iono_delay_klob = 0;
	}
	
	m_DelayCalc.iono_delay_klob = m_AtmDelayCalc.iono_delay_klob;
	m_DelayCalc.tropo_delay_hop = m_AtmDelayCalc.tropo_delay_hop;


}

/************************************************************************************************
* Function               : HopfieldModel
*
* Description            : This model computes the tropospheric delay of satellite signals.
*                          Reference: GPS Theory and Practice (3rd revised edition)
*                          by B.Hofmann-Wellenhof, H.Lichtenegger, and J.Collins
*
* Formal parameter(s)    : None
*
* Specific library calls : None
*
* Functions Called       : StdAtm()
*
* Version History:
* <1.1> <Thakur Shivam Singh> <10-10-2017>
*      1.1
****************************************************************************************/

void CDelayCalc::TropoDelayHop()
{
	double surface_dry_zenith_delay,
		surface_wet_zenith_delay;
	double max_dry_h;
	double alpha_dry[9],
		alpha_wet[9],
		r_wet,
		r_dry;

	double a_dry,
		b_dry,
		a_wet,
		b_wet,
		dry_sum,
		wet_sum,
		height_above_sea;

	double dry_delay,
		wet_delay;

	double const max_wet_h = 11000;
	double const c_one = 77.64e-2;
	double const c_two = -12.96e-2;
	double const c_three = 3.718e+3;

	// to find out
	double water_vap_pres,
		sat_vap_pressure;

	double some_val;

	height_above_sea = geo.ConvertHeight(m_Lat * 180 / EIGEN_PI, m_Longi * 180 / EIGEN_PI, m_Height, GeographicLib::Geoid::ELLIPSOIDTOGEOID);;

	// Call function to calcualte temperature, pressure and density using standard atmosphere model
	StdAtm(height_above_sea);

	// Calculate saturated and actual vapour pressure using
	// The ASCE Standardized Reference Evapotranspiration Equation
	some_val = (17.27 * (m_Temp - 273.15)) / ((m_Temp - 273.15) + 237.3);

	sat_vap_pressure =  610.8 * pow(Exp_Const,some_val);
	water_vap_pres =  (m_Humidity * sat_vap_pressure);

	// Calculate dry and wet zenith delay at the surface
	surface_dry_zenith_delay = c_one * (m_Pres / m_Temp);
	surface_wet_zenith_delay = (c_two * (water_vap_pres / m_Temp)) + (c_three * (water_vap_pres / pow(m_Temp, 2)));

	// Calculate height above user at which dry refractivity is zero
	max_dry_h = 40136 + 148.72*(m_Temp - 273.16);

	// Calculate r_dry and r_wet 6.110
	r_dry = sqrt(pow((RadiusOfEarth_Const + max_dry_h), 2) - pow((RadiusOfEarth_Const*cos(m_Ele)), 2)) -
		RadiusOfEarth_Const*sin(m_Ele);
	r_wet = sqrt(pow((RadiusOfEarth_Const + max_wet_h), 2) - pow((RadiusOfEarth_Const*cos(m_Ele)), 2)) -
		RadiusOfEarth_Const*sin(m_Ele);

	// calculate 'a' and 'b' parameters for dry and wet delays
	a_dry = -(sin(m_Ele) / max_dry_h);
	b_dry = -(pow(cos(m_Ele), 2) / (2 * RadiusOfEarth_Const*max_dry_h));

	a_wet = -(sin(m_Ele) / max_wet_h);
	b_wet = -(pow(cos(m_Ele), 2) / (2 * RadiusOfEarth_Const*max_wet_h));

	// Calculate alpha parameter for dry delay and wet delay
	alpha_dry[0] = 1;
	alpha_dry[1] = 4 * a_dry;
	alpha_dry[2] = (6 * pow(a_dry, 2)) + (4 * b_dry);
	alpha_dry[3] = 4 * a_dry*(pow(a_dry, 2) + (3 * b_dry));
	alpha_dry[4] = (pow(a_dry, 4)) + (12 * pow(a_dry, 2)*b_dry) + (6*pow(b_dry,2));
	alpha_dry[5] = 4 * a_dry*b_dry*(pow(a_dry, 2) + (3 * b_dry));
	alpha_dry[6] = pow(b_dry, 2)*((6 * pow(a_dry, 2)) + (4 * b_dry));
	alpha_dry[7] = 4 * a_dry*pow(b_dry, 3);
	alpha_dry[8] = pow(b_dry, 4);

	alpha_wet[0] = 1;
	alpha_wet[1] = 4 * a_wet;
	alpha_wet[2] = (6 * pow(a_wet, 2)) + (4 * b_wet);
	alpha_wet[3] = 4 * a_wet*(pow(a_wet, 2) + (3 * b_wet));
	alpha_wet[4] = (pow(a_wet, 4)) + (12 * pow(a_wet, 2)*b_wet) + (6 * pow(b_wet, 2));
	alpha_wet[5] = 4 * a_wet*b_wet*(pow(a_wet, 2) + (3 * b_wet));
	alpha_wet[6] = pow(b_wet, 2)*((6 * pow(a_wet, 2)) + (4 * b_wet));
	alpha_wet[7] = 4 * a_wet*pow(b_wet, 3);
	alpha_wet[8] = pow(b_wet, 4);

	wet_sum = 0;
	dry_sum = 0;
	for (int iter = 0; iter < 9; iter++)
	{
		dry_sum = dry_sum + ((alpha_dry[iter] / (iter + 1))*pow(r_dry, (iter + 1)));
		wet_sum = wet_sum + ((alpha_wet[iter] / (iter + 1))*pow(r_wet, (iter + 1)));

	}
	
	// Calculate dry and wet delays
	dry_delay = 1e-6*surface_dry_zenith_delay*dry_sum;
	wet_delay = 1e-6*surface_wet_zenith_delay*wet_sum;

	// Calculate tropospheric delay from modified Hopfield model
	m_AtmDelayCalc.tropo_delay_hop = dry_delay + wet_delay;
	

}

/************************************************************************************************
* Function               : StdAtm
*
* Description            : This model computes the Pressure, Temperature and Density
						   of the standard atmosphere at user position
*
* Formal parameter(s)    : None
*
* Return Value           : None
*
* Specific library calls : pow
*
* Functions Called       : None
*
* Assumption             : Maximum geopotential height of antenna is 30km
*
* Version History:
* <1.1> <Author> <THkaur Shivam Singh> <10/10/2017>
*      1.1
****************************************************************************************/

void CDelayCalc::StdAtm(double height_above_sea)
{
	// Geopotential Height
	double geop_height;

	// Temperature, pressure and density at different levels of atmosphere
	double t_one,
		t_two,
		p_one,
		p_two,
		rho_one,
		rho_two;
	
	// Calculate geopotential height from geometric height
	geop_height = (RadiusOfEarth_Const * height_above_sea) / (height_above_sea + RadiusOfEarth_Const);

	// Calculate temperature, pressure and density at height h1 of Standard Atmosphere model
	t_one = SeaLvlTemp_Const + (aNot_Const*hOne_Const);
	p_one = SeaLvlPre_Const*pow((t_one / SeaLvlTemp_Const), (-Gravity_Const / (aNot_Const*RSpecificAir_Const)));
	rho_one = SeaLvlDensity_Const*pow((t_one / SeaLvlTemp_Const), ((-Gravity_Const / (aNot_Const*RSpecificAir_Const)) - 1));

	// Calculate temperature, pressure and density at height h2 of Standard Atmosphere model
	t_two = t_one;
	p_two = p_one*pow(Exp_Const, -((Gravity_Const*(hTwo_Const - hOne_Const)) / (RSpecificAir_Const*t_two)));
	rho_two = rho_one*pow(Exp_Const, -((Gravity_Const*(hTwo_Const - hOne_Const)) / (RSpecificAir_Const*t_two)));

	// Calculate temperature, pressure and density at the antenna height
	if (geop_height < hOne_Const)
	{
		m_Temp = SeaLvlTemp_Const + (aNot_Const*geop_height);
		m_Pres = SeaLvlPre_Const*pow((m_Temp / SeaLvlTemp_Const), (-Gravity_Const / (aNot_Const*RSpecificAir_Const)));
		m_Density = SeaLvlDensity_Const*pow((m_Temp / SeaLvlTemp_Const), ((-Gravity_Const / (aNot_Const*RSpecificAir_Const)) - 1));
	}

	else if (geop_height <= hTwo_Const)
	{
		m_Temp = t_one;
		m_Pres = p_one * pow(Exp_Const, -((Gravity_Const*(geop_height - hOne_Const)) / (RSpecificAir_Const*m_Temp)));
		m_Density = rho_one * pow(Exp_Const, -((Gravity_Const*(geop_height - hOne_Const)) / (RSpecificAir_Const*m_Temp)));

	}

	else if (geop_height <= hThree_Const)
	{
		m_Temp = t_two + aTwo_Const*(geop_height - hTwo_Const);
		m_Pres = p_two * pow((m_Temp / t_two), (-Gravity_Const / (aTwo_Const*RSpecificAir_Const)));
		m_Density = rho_two * pow((m_Temp / SeaLvlTemp_Const), ((-Gravity_Const / (aTwo_Const*RSpecificAir_Const)) - 1));

	}

}

/************************************************************************************************
* Function               : IonoDelayKlob()
*
* Description            : This model computes the ionospheric delay of satellite signals using Klobuchar model.
*                          
*
* Formal parameter(s)    : None
*
* Return Value           : None
*
* System calls Used      : None
*
* Specific library calls : pow
*                          sin
*                          sqrt
*
* Functions Called       : None
*
* Version History        :
* <Version Number> <Author> <date> <defect Number> <Modification made and reason>
*      1.1
****************************************************************************************/

void CDelayCalc::IonoDelayKlob()
{
	double phi_i, psi, lambda_i, phi_m, solar_time,
			obliq_factor, period, amplitude, phase;
	double longi, azimuth, elevation, latitude;

	// alpha: coefficients of a cubic equation representing the amplitude 
	// of the vertical delay (4 coefficients - 8 bits each)
	// beta: coefficients of a cubic equation representing the period
	// of the model (4 coefficients - 8 bits each)
	//double alpha[4] = {0.1211E-07,  0.1490E-07, -0.5960E-07, -0.1192E-06};
	//double beta[4] = {0.9626E+05,  0.8192E+05, -0.1966E+06, -0.3932E+06};

	longi = m_Longi;
	azimuth = m_Azi;
	elevation = m_Ele;
	latitude = m_Lat;

	// Converting andle in radians to semi circles
	double semi2rad = Pi_Const;
	double rad2semi = 1 / semi2rad;

	// Convert azimuth to positive(East longitude)
	if (azimuth < 0)
	{
		azimuth = azimuth + 2 * Pi_Const;
	}

	// Find all elevation below 0.1 Deg and set them to 0.1 deg
	if (elevation < 0.0017)
	{
		elevation = 0.0017;
	}

	// Convert azimuth, elevation, latitude and longitude from radian to semicircles;
	elevation = elevation * rad2semi;
	azimuth = azimuth * rad2semi;
	latitude = latitude * rad2semi;
	longi = longi * rad2semi;

	// Calculate earth's central angle between the user position and the earth projection
	// of ionospheric intersection point (unit: semi-circles)
	psi = (0.0137 / (elevation + 0.11)) - 0.022;

	// Calculate geodetic latitude of the earth projection 
	// of the ionospheric intersection point (unit: semi-circles)
	phi_i = latitude + psi * cos(azimuth * semi2rad);
	if (phi_i > 0.416)
	{
		phi_i = 0.416;
	}
	if (phi_i < -0.416)
	{
		phi_i = -0.416;
	}

	// Calculate geodetic longitude of the earth projection 
	// of the ionospheric intersection point (unit: semi-circles)
	lambda_i = longi + ((psi * sin(azimuth * semi2rad)) / cos(phi_i * semi2rad));

	// Calculate geomagnetic latitude of the earth projection of the 
	// ionospheric intersection point (mean ionospheric height assumed 350 km) (unit: semi-circles)
	phi_m = phi_i + (0.064 * cos((lambda_i - 1.617) * semi2rad));

	// Calculate local solar time
	solar_time = (4.32 * (1E+4))* lambda_i + m_TimeOfWeek;
	solar_time = solar_time - int(solar_time / 86400) * 86400;

	if (solar_time > 86400)
	{
		solar_time = solar_time - 86400;
	}
	if (solar_time < 0)
	{
		solar_time = solar_time + 86400;
	}

	// Calculate obliquity factor
	obliq_factor = 1.0 +  (16.0 *pow((0.53 - elevation), 3));

	// Calculate the period
	period = m_Beta[0] + (m_Beta[1] * pow(phi_m, 1)) + (m_Beta[2] * pow(phi_m, 2)) + (m_Beta[3] * pow(phi_m, 3));
	if (period < 72000)
	{
		period = 72000;
	}

	// Calcuate the phase
	phase = (2 * Pi_Const*(solar_time - 50400)) / period;

	// Calcluate the amplitude 
	amplitude = m_Alpha[0] + (m_Alpha[1] * pow(phi_m, 1)) + (m_Alpha[2] * pow(phi_m, 2)) + 
					(m_Alpha[3] * pow(phi_m, 3));

	if (amplitude < 0)
	{
		amplitude = 0;
	}

	// Calcualte the ionospheric delay
	// daytime: -1.57 < phase < 1.57 
	// nightime: phase >= 1.57 or phase <= -1.57
	if (phase > -1.57 && phase < 1.57)
	{
		m_AtmDelayCalc.iono_delay_klob = (obliq_factor * ((5 * (1E-9)) + amplitude*(1 - (pow(phase, 2) / 2) +  \
										(pow(phase, 4) / 24))))*SpeedLight_Const;
	}
	else
	{
		m_AtmDelayCalc.iono_delay_klob = (obliq_factor * (5 * (1E-9)))*SpeedLight_Const;
	}

}