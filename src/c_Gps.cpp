/********************************************************************************************************
* Description            : Calculation of GPS Satellites ECEF position-velocity, range and range rates with respect to a
*`				           stationary user on Earth
* Specific library calls :  pow
*                           sqrt
*                           atan
*                           asin
* Classes                : 1. GPS
* Assumptions            : 1. rate of inclination of satellite orbit is zero
*                          2. Correction term for all sine and cosine harmonics are zero
* Reference              : 1. GPS Interface Control Document IS-GPS-200H
*                          2. Chapter 7, Paul D Grooves : Principles of GNSS, Inertial, and Multisensor Integrated Navigation
*                             System
* Version History        :
* <1.1><Thakur Shivam Singh><20/3/2017>
***********************************************************************************************************/

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<math.h>
#include<stdlib.h>
#include<Eigen/Dense>
using namespace std;
using namespace Eigen;
#include "c_Gps.h"
#include "ConstVar.h"
#include "c_main.h"



CGps::CGps(){

	// Allocate memory to struct pointer
	mp_Yuma = new GpsYumaData[MaxGpsSat_Const];
	memset(mp_Yuma, 0, sizeof(GpsYumaData)*MaxGpsSat_Const);

	mp_Gps = new Sat[MaxGpsSat_Const];
	memset(mp_Gps, 0, sizeof(Sat)*MaxGpsSat_Const);

	mp_GpsSat = new SatData[MaxGpsSat_Const];
	memset(mp_GpsSat, 0, sizeof(SatData)*MaxGpsSat_Const);

	mp_TransitTime = new double[MaxGpsSat_Const];
	memset(mp_TransitTime, 0, sizeof(double)*MaxGpsSat_Const);

}

CGps::~CGps(){
	
	// Deallocate memory to struct pointer
	delete(mp_Gps);
	mp_Gps = NULL;
	delete(mp_Yuma);
	mp_Yuma = NULL;
	delete(mp_GpsSat);
	mp_GpsSat = NULL;
	delete(mp_TransitTime);
	mp_TransitTime = NULL;
}

bool CGps::m_flag_gps = false;

void CGps::GpsSat(struct TimeVar m_TVar, struct UserMot m_UserMotion, struct SatData *mp_GpsSatData, double *mp_TransTime)
{

	
	// Allocate values from the function parameter
	memcpy(&m_GpsTVar, &m_TVar, sizeof(m_GpsTVar));
	memcpy(&m_GpsUserMotion, &m_UserMotion, sizeof(m_GpsUserMotion));
	memcpy(mp_TransitTime, mp_TransTime, sizeof(double)*MaxGpsSat_Const);


	/* function call to read YUMA file */
	while (!m_flag_gps){
		m_NoOfSat = YumaRead();
		m_flag_gps = true;
	}

	/* Calculate ECEF position and velocity of GPS satellites */
	for (int gps_sat_id = 0; gps_sat_id <= m_NoOfSat; gps_sat_id++)
	{
		// Calculation of time of interest  
		//m_TInt = m_GpsTVar.time_of_week + (m_GpsTVar.week_no - mp_Yuma[gps_sat_id].week)*SecsInWeek_Const ;

		// Calculate time difference between time of interest and time of applicability
		m_TimeDiff = m_GpsTVar.time_of_week + (m_GpsTVar.week_no - mp_Yuma[gps_sat_id].week)*SecsInWeek_Const
					- mp_Yuma[gps_sat_id].time;
		
		
		/* Calculate satellite position,velocity at signal receive time */
		PosiVel(gps_sat_id);

		/* Calculate range and range rate at signal receive time */
		RangeAndRate(gps_sat_id);

		// Recalculation is being done as we need the value of satellite parameters 
		// at the time of signal transmission
		// Calculation of transit time i.e time take by signal to reach from satellite to receiver
		mp_TransitTime[gps_sat_id] = mp_GpsSat[gps_sat_id].range / SpeedLight_Const;

		/* Caclulate signal transmission time */
		//m_SigTransmitTime = m_TInt - mp_TransitTime[gps_sat_id];

		// Calculate time differnece as the difference between time of applicability and time of signal transmission
		m_TimeDiff = m_GpsTVar.time_of_week + (m_GpsTVar.week_no - mp_Yuma[gps_sat_id].week)*SecsInWeek_Const
					- mp_Yuma[gps_sat_id].time - mp_TransitTime[gps_sat_id];

		/*
		// Update time difference so as to account for week crossover
		if (m_time_diff > 302400)
		{
			m_time_diff -= 604800;
		}
		if (m_time_diff < -302400)
		{
			m_time_diff += 604800;
		}
		*/
		// Caclculate satellite position,velocity at signal transmission time in ecef at recieve time
		PosiVel(gps_sat_id);
		// Calculate range and range rate at signal transmission time 
		RangeAndRate(gps_sat_id);

		// Caclculate azimuth and elevation angle between user and satellite at signal transmission time 
		AziEle(gps_sat_id);


		// Apply clock correction term if satellite is visible i.e elevation > mask angle
		// mask angle is given by the user, currently putting its value manually 
		if (mp_GpsSat[gps_sat_id].elevation > 0.0872665)
		{
			mp_GpsSat[gps_sat_id].clock_correction = SpeedLight_Const*(mp_Yuma[gps_sat_id].af0s + mp_Yuma[gps_sat_id].af1s *
				(int(m_GpsTVar.time_of_day) % 7200));
			mp_GpsSat[gps_sat_id].visible = true;
		}
		else
		{
			mp_GpsSat[gps_sat_id].clock_correction = 0;
			mp_GpsSat[gps_sat_id].visible = false;
		}

		/*
		// Check whether the elevation is greater than threshold or not
		if (mp_GPS_sat[gps_sat_id].elevation > Limit_Ele_Const)
		{
			mp_gps_sat_data[gps_sat_id].elevation_check = 1;
		}
		*/
			
	}

	// Pass the value of calculated GPS satellite variables back to the main class
	memcpy(mp_GpsSatData, mp_GpsSat, sizeof(SatData)*MaxGpsSat_Const);
	memcpy(mp_TransTime, mp_TransitTime, sizeof(double)*MaxGpsSat_Const);

}


/********************************************************************************************************
* Function               : YumaRead
* Description            : Function to read Yuma file
* Specific library calls :  None
* Assumptions            :  None
* Reference              :  None
* Version History        :
* <1.1><Thakur Shivam Singh><20/3/2017>
***********************************************************************************************************/

int CGps::YumaRead()
{

	// take input from file 
	ifstream iFile("gps_yuma.txt");

	// Output file 
	ofstream ofile_1("gps_orbital.txt");

	string strn;
	double inte;

	// Commands to read the Yuma file, reading word by word in each line 
	int count = 0;
	int index = -1;
	while (!iFile.eof())
	{
		string temp;
		iFile >> temp;
		if (iFile.eof()) break;

		if (temp == "********")
		{
			iFile >> strn; iFile >> strn; iFile >> strn; iFile >> strn; iFile >> strn; iFile >> strn;
			index++;
			count = 0;
			continue;
		}
		if (temp == "ID:")
		{
			iFile >> inte;
			mp_Yuma[index].sat_id = inte;
		}
		else if (temp == "Health:")
		{
			iFile >> inte;
			mp_Yuma[index].health = inte;
		}
		else if (temp == "Eccentricity:")
		{
			iFile >> inte;
			mp_Yuma[index].eccent = inte;
		}
		else if (temp == "Time")
		{
			iFile >> strn; iFile >> strn; iFile >> inte;
			mp_Yuma[index].time = inte;
		}
		else if (temp == "Orbital")
		{
			iFile >> strn; iFile >> inte;
			mp_Yuma[index].a_o_i = inte;
		}
		else if (temp == "Rate")
		{
			iFile >> strn; iFile >> strn; iFile >> strn; iFile >> inte;
			mp_Yuma[index].rt_r_asc = inte;
		}
		else if (temp == "SQRT(A)")
		{
			iFile >> strn; iFile >> strn; iFile >> inte;

			mp_Yuma[index].sqrt_semi_maj = inte;
			mp_Yuma[index].semi_maj = (inte*inte);
		}
		else if (temp == "Right")
		{
			iFile >> strn; iFile >> strn; iFile >> strn; iFile >> inte;
			mp_Yuma[index].r_asc = inte;
		}
		else if (temp == "Argument")
		{
			iFile >> strn; iFile >> strn; iFile >> inte;
			mp_Yuma[index].arg_of_pge = inte;
		}
		else if (temp == "Mean")
		{
			iFile >> strn;
			iFile >> inte;
			mp_Yuma[index].m_ini = inte;
		}
		else if (temp == "Af0(s):")
		{
			iFile >> inte;
			mp_Yuma[index].af0s = inte;
		}
		else if (temp == "Af1(s/s):")
		{

			iFile >> inte;
			mp_Yuma[index].af1s = inte;
		}
		else if (temp == "week:")
		{

			iFile >> inte;
			mp_Yuma[index].week = inte;
		}

	}
	// Save the yuma file data in a file 
	for (int iter = 0; iter <= index; iter++)
	{
		ofile_1 << mp_Yuma[iter].sat_id << "  " << mp_Yuma[iter].health << "  " << mp_Yuma[iter].eccent
			<< "  " << mp_Yuma[iter].time << "  " << mp_Yuma[iter].a_o_i << "  " << mp_Yuma[iter].rt_r_asc
			<< "  " << mp_Yuma[iter].semi_maj << "  " << mp_Yuma[iter].r_asc << "  " << mp_Yuma[iter].arg_of_pge
			<< "  " << mp_Yuma[iter].m_ini << "  " << mp_Yuma[iter].af0s << "  " << mp_Yuma[iter].af1s << "  "
			<< mp_Yuma[iter].week << endl;
		ofile_1 << "#$$$$$\n";

	}
	return index;
}
/*******************************************************************************************************/

/********************************************************************************************************
* Function               : PosiVel
* Description            : Calculation of GPS Satellites ECEF position-velocity,
* Specific library calls :  pow
*                           sqrt
*                           atan
*                           asin
* Function parameter     : gps_sat_id
* Return value           : none
* Assumptions            : 1. rate of inclination of satellite orbit is zero
*                          2. Correction term for all sine and cosine harmonics are zero
* Reference              : 1. GPS Interface Control Document IS-GPS-200H
*                          2. Chapter 7, Paul D Grooves : Principles of GNSS, Inertial, and Multisensor Integrated Navigation
*                              System
* Version History        :
* <1.1><Thakur Shivam Singh><20/3/2017>
***********************************************************************************************************/


void CGps::PosiVel(int gps_sat_id)
{
	double t_amoly, e_amoly, e_one, eccent_amoly;
	double mean_motion;
	double x_pos, y_pos, x_pos_dot, y_pos_dot;
	double c_m_amoly, arg_of_lat, c_rad, c_longi_asc_node;
	double eccent_amoly_dot, arg_of_lat_dot, c_rad_dot, c_longi_asc_node_dot;


	// Mean Motion 
	mean_motion = sqrt(GMProd_Const / pow(mp_Yuma[gps_sat_id].semi_maj, 3));

	// Corrected Mean Anomaly  
	c_m_amoly = mp_Yuma[gps_sat_id].m_ini + mean_motion*m_TimeDiff;

	// Eccentric Anomaly from iterative process 
	e_one = c_m_amoly;
	do
	{
		eccent_amoly = e_one;
		e_one = eccent_amoly - ((eccent_amoly - mp_Yuma[gps_sat_id].eccent*sin(eccent_amoly)
			- c_m_amoly) / (1 - mp_Yuma[gps_sat_id].eccent*cos(eccent_amoly)));
	} while ((e_one - eccent_amoly) > 0.0001);

	e_amoly = e_one;

	// true anomaly from eccentric anomaly 
	t_amoly = atan2(((sqrt(1 - pow(mp_Yuma[gps_sat_id].eccent, 2)))*sin(e_amoly)),
		(cos(e_amoly) - mp_Yuma[gps_sat_id].eccent));

	// argument of latitude 
	arg_of_lat = t_amoly + mp_Yuma[gps_sat_id].arg_of_pge;

	// Corrected radius 
	c_rad = mp_Yuma[gps_sat_id].semi_maj*(1 - mp_Yuma[gps_sat_id].eccent*cos(e_amoly));

	// Position in orbital plane 
	x_pos = c_rad*cos(arg_of_lat);
	y_pos = c_rad*sin(arg_of_lat);

	
	// Corrected longitude of ascending node
	// Correction is given so as to properly account for the ECEF frame 
	// at time of signal transmission and receive time
	c_longi_asc_node = mp_Yuma[gps_sat_id].r_asc + mp_Yuma[gps_sat_id].rt_r_asc*m_TimeDiff -
		RtEarthRotn_Const*(m_TimeDiff + mp_Yuma[gps_sat_id].time + mp_TransitTime[gps_sat_id]);

	// Rate of change of different variables 
	/**********************************************************************************************************/
	// Rate of change of eccentric anomaly Equation : 7.16 
	eccent_amoly_dot = mean_motion / (1 - mp_Yuma[gps_sat_id].eccent*cos(e_amoly));

	// Rate of change of argument of latitude Equation : 7.17 
	arg_of_lat_dot = (sin(t_amoly)* eccent_amoly_dot) / sin(e_amoly);

	// Rate of change of orbital radius Equation : 7.18 
	c_rad_dot = mp_Yuma[gps_sat_id].semi_maj*mp_Yuma[gps_sat_id].eccent*sin(e_amoly)*eccent_amoly_dot;

	// Rate of change of position in orbital plane Equation : 7.19 
	x_pos_dot = c_rad_dot*cos(arg_of_lat)
		- c_rad*arg_of_lat_dot*sin(arg_of_lat);


	y_pos_dot = c_rad_dot*sin(arg_of_lat)
		+ c_rad*arg_of_lat_dot*cos(arg_of_lat);


	// Rate of change of longitude of ascending node Equation : 7.20 
	c_longi_asc_node_dot = mp_Yuma[gps_sat_id].rt_r_asc - RtEarthRotn_Const;

	/****************************************************************************************************************/

	// ECEF Position-velocity of satellites  
	// Earth fixed coordinates
	mp_GpsSat[gps_sat_id].x_cord = x_pos*cos(c_longi_asc_node)
		- y_pos*cos(mp_Yuma[gps_sat_id].a_o_i)*sin(c_longi_asc_node);

	mp_GpsSat[gps_sat_id].y_cord = x_pos*sin(c_longi_asc_node)
		+ y_pos*cos(mp_Yuma[gps_sat_id].a_o_i)*cos(c_longi_asc_node);

	mp_GpsSat[gps_sat_id].z_cord = y_pos*sin(mp_Yuma[gps_sat_id].a_o_i);

	// ECEF Velocity of satellites Equation : 7.22  
	mp_GpsSat[gps_sat_id].x_vel = (x_pos_dot* cos(c_longi_asc_node)
		- y_pos_dot*cos(mp_Yuma[gps_sat_id].a_o_i)*sin(c_longi_asc_node))
		- (c_longi_asc_node_dot*(x_pos* sin(c_longi_asc_node)
		+ y_pos*cos(mp_Yuma[gps_sat_id].a_o_i)*cos(c_longi_asc_node)));

	mp_GpsSat[gps_sat_id].y_vel = (x_pos_dot* sin(c_longi_asc_node)
		+ y_pos_dot*cos(mp_Yuma[gps_sat_id].a_o_i)*cos(c_longi_asc_node))
		- (c_longi_asc_node_dot*(-x_pos* cos(c_longi_asc_node)
		+ y_pos*cos(mp_Yuma[gps_sat_id].a_o_i)*sin(c_longi_asc_node)));

	mp_GpsSat[gps_sat_id].z_vel = y_pos_dot * sin(mp_Yuma[gps_sat_id].a_o_i);
}

/*****************************************************************************************************************/


/********************************************************************************************************
* Function               : RangeAndRate
* Description            : Calculation of range and range rates with respect to a
*`				           stationary user on Earth
* Specific library calls :  pow
*                           sqrt
*                           atan
*                           asin
* Function Parameter     : gps_sat_id
* Return                 : none
* Assumptions            : 1. rate of inclination of satellite orbit is zero
*                          2. Correction term for all sine and cosine harmonics are zero
* Reference              :  1. GPS Interface Control Document IS-GPS-200H
*                           2. Chapter 7, Paul D Grooves : Principles of GNSS, Inertial,
*                              and Multisensor Integrated Navigation Systems
* Version History        :
* <1.1><Thakur Shivam Singh><20/3/2017>
***********************************************************************************************************/


void CGps::RangeAndRate(int gps_sat_id)
{

	double temp_one;
	double earth_rotn_corr, earth_rotn_corr_rt;
	// Calculate of Range and Range rate 
	// Earth-rotation correction Equation : 7.27 
	earth_rotn_corr = (RtEarthRotn_Const / SpeedLight_Const)*(mp_GpsSat[gps_sat_id].y_cord*(m_GpsUserMotion.user_pos[0])
						- mp_GpsSat[gps_sat_id].x_cord*(m_GpsUserMotion.user_pos[1]));

	// Range calculation  Equation : 7.26
	temp_one = sqrt(pow((mp_GpsSat[gps_sat_id].x_cord - m_GpsUserMotion.user_pos[0]), 2)
					+ pow((mp_GpsSat[gps_sat_id].y_cord - m_GpsUserMotion.user_pos[1]), 2)
					+ pow((mp_GpsSat[gps_sat_id].z_cord - (m_GpsUserMotion.user_pos[2])), 2));

	mp_GpsSat[gps_sat_id].range = temp_one;

	// Rate of Earth-rotation correction Equation : 7.40  
	earth_rotn_corr_rt = (RtEarthRotn_Const / SpeedLight_Const)*(mp_GpsSat[gps_sat_id].y_vel*m_GpsUserMotion.user_pos[0]
						- mp_GpsSat[gps_sat_id].x_vel*m_GpsUserMotion.user_pos[1]);

	// Line of sight vector from user to satellite Equation 7.35
	mp_Gps[gps_sat_id].x_unit = (mp_GpsSat[gps_sat_id].x_cord - m_GpsUserMotion.user_pos[0]) / (temp_one) ;
	mp_Gps[gps_sat_id].y_unit = (mp_GpsSat[gps_sat_id].y_cord - m_GpsUserMotion.user_pos[1]) / (temp_one);
	mp_Gps[gps_sat_id].z_unit = (mp_GpsSat[gps_sat_id].z_cord - m_GpsUserMotion.user_pos[2]) / (temp_one);


	// Calculation of range rate Equation : 7.39 
	mp_GpsSat[gps_sat_id].range_rt = mp_Gps[gps_sat_id].x_unit*(mp_GpsSat[gps_sat_id].x_vel - m_GpsUserMotion.user_vel[0])
									+ mp_Gps[gps_sat_id].y_unit*(mp_GpsSat[gps_sat_id].y_vel - m_GpsUserMotion.user_vel[1])
									+ mp_Gps[gps_sat_id].z_unit*(mp_GpsSat[gps_sat_id].z_vel - m_GpsUserMotion.user_vel[2]);

}

/********************************************************************************************************
* Function               : AziEle
*
* Description            : Function to calculate Azimuth and Elevation angles
*
* Function Parameter     : gps_sat_id
*
* Return value           : None
*
* Specific library calls : None
*
* Functions called       : None
*
* Reference              : 1. GPS Interface Control Document IS-GPS-200H
*						   2. Chapter 7, Paul D Grooves : Principles of GNSS, Inertial, and Multisensor Integrated Navigation
*						      Systems
*
*Version History        :
* <1.1><Thakur Shivam Singh><20/3/2017>
***********************************************************************************************************/

void CGps::AziEle(int gps_sat_id)
{

	Matrix3d tran_ecef_local(3, 3);
	MatrixXd lineofsight_ecef(3, 1);
	MatrixXd lineofsight_local(3, 1);
	
	// Line of sight in ECEF frame 
	lineofsight_ecef << mp_Gps[gps_sat_id].x_unit,
		mp_Gps[gps_sat_id].y_unit,
		mp_Gps[gps_sat_id].z_unit;

	// Coordinate transformation matrix from ECEF to Local Navigation Frame
	tran_ecef_local << -sin(m_GpsUserMotion.lat)*cos(m_GpsUserMotion.longi), \
						-sin(m_GpsUserMotion.lat)*sin(m_GpsUserMotion.longi), cos(m_GpsUserMotion.lat),
					-sin(m_GpsUserMotion.longi), cos(m_GpsUserMotion.longi), 0,
					-cos(m_GpsUserMotion.lat)*cos(m_GpsUserMotion.longi), \
					-cos(m_GpsUserMotion.lat)*sin(m_GpsUserMotion.longi), -sin(m_GpsUserMotion.lat);
	/***************************************************************************/

	// Transform line of sight vector in ECEf frame to local navigation frame 
	lineofsight_local = tran_ecef_local*lineofsight_ecef;

	// Calculate azimuth and elevation angles 
	mp_GpsSat[gps_sat_id].elevation = -(asin(lineofsight_local(2, 0)));
	mp_GpsSat[gps_sat_id].azimuth = (atan2(lineofsight_local(1, 0), lineofsight_local(0, 0)));
}
