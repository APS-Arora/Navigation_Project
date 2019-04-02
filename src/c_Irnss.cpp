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
#include "c_Irnss.h"
#include "ConstVar.h"
using namespace std;
using namespace Eigen;
#include "c_main.h"


CIrnss::CIrnss(){

	// Allocate memory to struct pointer
	mp_IrnssYuma = new IrnssYumaData[MaxIrnssSat_Const];
	memset(mp_IrnssYuma, 0, sizeof(IrnssYumaData)*MaxIrnssSat_Const);

	mp_Irnss = new Sat[MaxIrnssSat_Const];
	memset(mp_Irnss, 0, sizeof(Sat)*MaxIrnssSat_Const);

	mp_IrnssSat = new SatData[MaxIrnssSat_Const];
	memset(mp_IrnssSat, 0, sizeof(SatData)*MaxIrnssSat_Const);

	mp_TransitTime = new double[MaxIrnssSat_Const];
	memset(mp_TransitTime, 0, sizeof(double)*MaxIrnssSat_Const);

}

CIrnss::~CIrnss(){

	// Deallocate memory to struct pointer
	delete(mp_Irnss);
	mp_Irnss = NULL;
	delete(mp_IrnssYuma);
	mp_IrnssYuma = NULL;
	delete(mp_IrnssSat);
	mp_IrnssSat = NULL;
	delete(mp_TransitTime);
	mp_TransitTime = NULL;
}

bool CIrnss::m_flag_irnss = false;

void CIrnss::IrnssSat(struct TimeVar m_TVar, struct UserMot m_UserMotion, struct SatData *mp_IrnssSatData, double *mp_IrnssTransTime)
{

	
	// Allocate values from the function parameter
	memcpy(&m_IrnssTVar, &m_TVar, sizeof(m_IrnssTVar));
	memcpy(&m_IrnssUserMotion, &m_UserMotion, sizeof(m_IrnssUserMotion));
	memcpy(mp_TransitTime, mp_IrnssTransTime, sizeof(double)*MaxIrnssSat_Const);


	/* function call to read YUMA file */
	while (!m_flag_irnss)
	{
		m_NoOfIrnssSat = IrnssYumaRead();
		m_flag_irnss = true;
	}

	/* calculation of ECEF position and velocity of GPS satellites */
	for (int irnss_sat_id = 0; irnss_sat_id < MaxIrnssSat_Const; irnss_sat_id++)
	{

		// Calculate time difference between time of interest and time of applicability
		m_TimeDiff = m_IrnssTVar.time_of_week + (m_IrnssTVar.week_no - mp_IrnssYuma[irnss_sat_id].week)*SecsInWeek_Const
			- mp_IrnssYuma[irnss_sat_id].time;

		/* Calculate satellite position,velocity at signal receive time */
		IrnssPosiVel(irnss_sat_id);

		/* Calculate range and range rate at signal receive time */
		IrnssRangeAndRate(irnss_sat_id);

		// Recalculation if being done as we need the value of satellite parameters 
		// at the time of signal transmission
		// Calculation of transit time i.e time take by signal to reach from satellite to receiver
		mp_TransitTime[irnss_sat_id] = (mp_IrnssSat[irnss_sat_id].range / SpeedLight_Const);

		// Calculate time difference between time of interest and time of applicability
		m_TimeDiff = m_IrnssTVar.time_of_week + (m_IrnssTVar.week_no - mp_IrnssYuma[irnss_sat_id].week)*SecsInWeek_Const
			- mp_IrnssYuma[irnss_sat_id].time - mp_TransitTime[irnss_sat_id];

		/* Caclculate satellite position,velocity at the time of signal transmission*/
		IrnssPosiVel(irnss_sat_id);

		// Calculate range and range rate at signal transmission time 
		IrnssRangeAndRate(irnss_sat_id);

		// Caclculate azimuth and elevation angle between user and satellite at signal transmission time
		IrnssAziEle(irnss_sat_id);
		
		// Apply clock correction term
		mp_IrnssSat[irnss_sat_id].clock_correction = SpeedLight_Const*(mp_IrnssYuma[irnss_sat_id].af0s + mp_IrnssYuma[irnss_sat_id].af1s *
															(int(m_IrnssTVar.time_of_day) % 7200));


		/*
		// Check whether the elevation angle is greater than the threshhold or not
		if (mp_IRNSS_sat[irnss_sat_id].elevation > Limit_Ele_Const)
		{
			mp_irnss_sat_data[irnss_sat_id][iter_no].elevation_check = 1;
		}
		*/
		
	}

	// Pass the value of calculated GPS satellite variables back to the main class
	memcpy(mp_IrnssSatData, mp_IrnssSat, sizeof(SatData)*MaxIrnssSat_Const);
	memcpy(mp_IrnssTransTime, mp_TransitTime, sizeof(double)*MaxIrnssSat_Const);

}


/********************************************************************************************************
* Function               : IrnssYumaRead
* Description            : Function to read Yuma file
* Specific library calls :  None
* Assumptions            :  None
* Reference              :  None
* Version History        :
* <1.1><Thakur Shivam Singh><20/3/2017>
***********************************************************************************************************/

int CIrnss::IrnssYumaRead()
{

	/* taking input from file */
	ifstream iFile("irnss_yuma.txt");

	/* Output file */
	ofstream ofile_1("Irnss_orbital.txt");

	string strn;
	double inte;

	/* Commands to read the Yuma file, reading word by word in each line */
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
			mp_IrnssYuma[index].sat_id = inte;
		}
		else if (temp == "Health:")
		{
			iFile >> inte;
			mp_IrnssYuma[index].health = inte;
		}
		else if (temp == "Eccentricity:")
		{
			iFile >> inte;
			mp_IrnssYuma[index].eccent = inte;
		}
		else if (temp == "Time")
		{
			iFile >> strn; iFile >> strn; iFile >> inte;
			mp_IrnssYuma[index].time = inte;
		}
		else if (temp == "Orbital")
		{
			iFile >> strn; iFile >> inte;
			mp_IrnssYuma[index].a_o_i = inte;
		}
		else if (temp == "Rate")
		{
			iFile >> strn; iFile >> strn; iFile >> strn; iFile >> inte;
			mp_IrnssYuma[index].rt_r_asc = inte;
		}
		else if (temp == "SQRT(A)")
		{
			iFile >> strn; iFile >> strn; iFile >> inte;

			mp_IrnssYuma[index].sqrt_semi_maj = inte;
			mp_IrnssYuma[index].semi_maj = (inte*inte);
		}
		else if (temp == "Right")
		{
			iFile >> strn; iFile >> strn; iFile >> strn; iFile >> inte;
			mp_IrnssYuma[index].r_asc = inte;
		}
		else if (temp == "Argument")
		{
			iFile >> strn; iFile >> strn; iFile >> inte;
			mp_IrnssYuma[index].arg_of_pge = inte;
		}
		else if (temp == "Mean")
		{
			iFile >> strn;
			iFile >> inte;
			mp_IrnssYuma[index].m_ini = inte;
		}
		else if (temp == "Af0(s):")
		{
			iFile >> inte;
			mp_IrnssYuma[index].af0s = inte;
		}
		else if (temp == "Af1(s/s):")
		{

			iFile >> inte;
			mp_IrnssYuma[index].af1s = inte;
		}
		else if (temp == "week:")
		{

			iFile >> inte;
			mp_IrnssYuma[index].week = inte;
		}

	}
	/* Saving the yuma file data in a file */
	for (int iter = 0; iter <= index; iter++)
	{
		ofile_1 << mp_IrnssYuma[iter].sat_id << "  " << mp_IrnssYuma[iter].health << "  " 
			<< mp_IrnssYuma[iter].eccent << "  " << mp_IrnssYuma[iter].time << "  " 
			<< mp_IrnssYuma[iter].a_o_i << "  " << mp_IrnssYuma[iter].rt_r_asc << "  " 
			<< mp_IrnssYuma[iter].semi_maj << "  " << mp_IrnssYuma[iter].r_asc << "  " 
			<< mp_IrnssYuma[iter].arg_of_pge << "  " << mp_IrnssYuma[iter].m_ini << "  " 
			<< mp_IrnssYuma[iter].af0s << "  " << mp_IrnssYuma[iter].af1s << "  " << mp_IrnssYuma[iter].week << endl;
		ofile_1 << "#$$$$$\n";

	}
	return index;
}
/*******************************************************************************************************/

/********************************************************************************************************
* Function               : IrnssPosiVel
* Description            : Calculation of IRNSS Satellites ECEF position-velocity,
* Specific library calls :  pow
*                           sqrt
*                           atan
*                           asin
* Function parameter     : iter
* Return value           : none
* Assumptions            : 1. rate of inclination of satellite orbit is zero
*                          2. Correction term for all sine and cosine harmonics are zero
* Reference              : 1. GPS Interface Control Document IS-GPS-200H
*                          2. Chapter 7, Paul D Grooves : Principles of GNSS, Inertial, and Multisensor Integrated Navigation
*                              System
* Version History        :
* <1.1><Thakur Shivam Singh><20/3/2017>
***********************************************************************************************************/


void CIrnss::IrnssPosiVel(int sat_id)
{
	double t_amoly, e_amoly, e_one, eccent_amoly;
	double mean_motion;
	double x_pos, y_pos, x_pos_dot, y_pos_dot;
	double c_m_amoly, arg_of_lat, c_rad, c_longi_asc_node;
	double eccent_amoly_dot, arg_of_lat_dot, c_rad_dot, c_longi_asc_node_dot;


	/* Mean Motion */
	mean_motion = sqrt(GMProd_Const / pow(mp_IrnssYuma[sat_id].semi_maj, 3));

	/* Corrected Mean Anomaly  */
	c_m_amoly = mp_IrnssYuma[sat_id].m_ini + mean_motion*m_TimeDiff;

	/* Eccentric Anomaly from iterative process */
	e_one = c_m_amoly;
	do
	{
		eccent_amoly = e_one;
		e_one = eccent_amoly - ((eccent_amoly - mp_IrnssYuma[sat_id].eccent*sin(eccent_amoly)
			- c_m_amoly) / (1 - mp_IrnssYuma[sat_id].eccent*cos(eccent_amoly)));
	} while ((e_one - eccent_amoly) > 0.0001);

	e_amoly = e_one;

	// true anomaly from eccentric anomaly 
	t_amoly = atan2(((sqrt(1 - pow(mp_IrnssYuma[sat_id].eccent, 2)))*sin(e_amoly)),
		(cos(e_amoly) - mp_IrnssYuma[sat_id].eccent));

	/*  Argument of latitude */
	arg_of_lat = t_amoly + mp_IrnssYuma[sat_id].arg_of_pge;

	/* Corrected radius */
	c_rad = mp_IrnssYuma[sat_id].semi_maj*(1 - mp_IrnssYuma[sat_id].eccent*cos(e_amoly));

	/* Position in orbital plane */
	x_pos = c_rad*cos(arg_of_lat);
	y_pos = c_rad*sin(arg_of_lat);


	/* Corrected longitude of ascending node */
	c_longi_asc_node = mp_IrnssYuma[sat_id].r_asc + mp_IrnssYuma[sat_id].rt_r_asc*m_TimeDiff -
		RtEarthRotn_Const*(m_TimeDiff + mp_TransitTime[sat_id] + mp_IrnssYuma[sat_id].time);

	/* Rate of change of different variables */
	/**********************************************************************************************************/
	/*  Rate of change of eccentric anomaly Equation : 7.16 */
	eccent_amoly_dot = mean_motion / (1 - mp_IrnssYuma[sat_id].eccent*cos(e_amoly));

	/* Rate of change of argument of latitude Equation : 7.17 */
	arg_of_lat_dot = (sin(t_amoly)* eccent_amoly_dot) / sin(e_amoly);

	/* Rate of change of orbital radius Equation : 7.18 */
	c_rad_dot = mp_IrnssYuma[sat_id].semi_maj*mp_IrnssYuma[sat_id].eccent*sin(e_amoly)*eccent_amoly_dot;

	/* Rate of change of position in orbital plane Equation : 7.19 */
	x_pos_dot = c_rad_dot*cos(arg_of_lat)
		- c_rad*arg_of_lat_dot*sin(arg_of_lat);


	y_pos_dot = c_rad_dot*sin(arg_of_lat)
		+ c_rad*arg_of_lat_dot*cos(arg_of_lat);

	/* Rate of change of longitude of ascending node Equation : 7.20 */
	c_longi_asc_node_dot = mp_IrnssYuma[sat_id].rt_r_asc - RtEarthRotn_Const;

	/****************************************************************************************************************/

	/* ECEF Position-velocity of satellites  */
	/* Earth fixed coordinates*/
	mp_IrnssSat[sat_id].x_cord = x_pos*cos(c_longi_asc_node)
		- y_pos*cos(mp_IrnssYuma[sat_id].a_o_i)*sin(c_longi_asc_node);

	mp_IrnssSat[sat_id].y_cord = x_pos*sin(c_longi_asc_node)
		+ y_pos*cos(mp_IrnssYuma[sat_id].a_o_i)*cos(c_longi_asc_node);

	mp_IrnssSat[sat_id].z_cord = y_pos*sin(mp_IrnssYuma[sat_id].a_o_i);

	/* ECEF Velocity of satellites Equation : 7.22  */
	mp_IrnssSat[sat_id].x_vel = (x_pos_dot* cos(c_longi_asc_node)
		- y_pos_dot*cos(mp_IrnssYuma[sat_id].a_o_i)*sin(c_longi_asc_node))
		- (c_longi_asc_node_dot*(x_pos* sin(c_longi_asc_node)
		+ y_pos*cos(mp_IrnssYuma[sat_id].a_o_i)*cos(c_longi_asc_node)));

	mp_IrnssSat[sat_id].y_vel = (x_pos_dot* sin(c_longi_asc_node)
		+ y_pos_dot*cos(mp_IrnssYuma[sat_id].a_o_i)*cos(c_longi_asc_node))
		- (c_longi_asc_node_dot*(-x_pos* cos(c_longi_asc_node)
		+ y_pos*cos(mp_IrnssYuma[sat_id].a_o_i)*sin(c_longi_asc_node)));

	mp_IrnssSat[sat_id].z_vel = y_pos_dot * sin(mp_IrnssYuma[sat_id].a_o_i);

}

/*****************************************************************************************************************/


/********************************************************************************************************
* Function               : IrnssRangeAndRate
* Description            : Calculation of range and range rates with respect to a
*`				           stationary user on Earth
* Specific library calls :  pow
*                           sqrt
*                           atan
*                           asin
* Function Parameter     : iter
* Return                 : none
* Assumptions            : 1. rate of inclination of satellite orbit is zero
*                          2. Correction term for all sine and cosine harmonics are zero
* Reference              : 1. Chapter 7, Paul D Grooves : Principles of GNSS, Inertial, and Multisensor Integrated Navigation
*                              Systems                        
* Version History        :
* <1.1><Thakur Shivam Singh><20/3/2017>
***********************************************************************************************************/


void CIrnss::IrnssRangeAndRate(int sat_id)
{
	double temp_one;

	/* Calculate of Range and Range rate */
	/* Range calculation  Equation : 7.26*/
	temp_one = sqrt(pow((mp_IrnssSat[sat_id].x_cord - m_IrnssUserMotion.user_pos[0]), 2)
				+ pow((mp_IrnssSat[sat_id].y_cord - m_IrnssUserMotion.user_pos[1]), 2)
				+ pow((mp_IrnssSat[sat_id].z_cord - m_IrnssUserMotion.user_pos[2]), 2));

	mp_IrnssSat[sat_id].range = temp_one ;

	/* Line of sight vector from user to satellite Equation 7.35*/
	mp_Irnss[sat_id].x_unit = (mp_IrnssSat[sat_id].x_cord - m_IrnssUserMotion.user_pos[0]) / temp_one;
	mp_Irnss[sat_id].y_unit = (mp_IrnssSat[sat_id].y_cord - m_IrnssUserMotion.user_pos[1]) / temp_one;
	mp_Irnss[sat_id].z_unit = (mp_IrnssSat[sat_id].z_cord - m_IrnssUserMotion.user_pos[2]) / temp_one;


	/* Calculation of range rate Equation : 7.39 */
	mp_IrnssSat[sat_id].range_rt = mp_Irnss[sat_id].x_unit*(mp_IrnssSat[sat_id].x_vel - m_IrnssUserMotion.user_vel[0])
									+ mp_Irnss[sat_id].y_unit*(mp_IrnssSat[sat_id].y_vel - m_IrnssUserMotion.user_vel[1])
									+ mp_Irnss[sat_id].z_unit*(mp_IrnssSat[sat_id].z_vel - m_IrnssUserMotion.user_vel[2]);

}

/********************************************************************************************************
* Function               : IrnssAziEle
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
* Reference              : 1. Chapter 7, Paul D Grooves : Principles of GNSS, Inertial,
*						      and Multisensor Integrated Navigation Systems			   
* Version History        :
* <1.1><Thakur Shivam Singh><20/3/2017>
***********************************************************************************************************/



void CIrnss::IrnssAziEle(int sat_id)
{
	Matrix3d tran_ecef_local(3, 3);
	MatrixXd lineofsight_ecef(3, 1);
	MatrixXd lineofsight_local(3, 1);

	/* Line of sight in ECEF frame */
	lineofsight_ecef << mp_Irnss[sat_id].x_unit,
		mp_Irnss[sat_id].y_unit,
		mp_Irnss[sat_id].z_unit;

	/* Coordinate transformation matrix from ECEF to Local Navigation Frame*/
	tran_ecef_local << -sin(m_IrnssUserMotion.lat)*cos(m_IrnssUserMotion.longi), \
					-sin(m_IrnssUserMotion.lat)*sin(m_IrnssUserMotion.longi), cos(m_IrnssUserMotion.lat),
				-sin(m_IrnssUserMotion.longi), cos(m_IrnssUserMotion.longi), 0,
				-cos(m_IrnssUserMotion.lat)*cos(m_IrnssUserMotion.longi), \
				-cos(m_IrnssUserMotion.lat)*sin(m_IrnssUserMotion.longi), -sin(m_IrnssUserMotion.lat);
	/***************************************************************************/

	/* Transformation of line of sight vector in ECEf frame to local navigation frame */
	lineofsight_local = tran_ecef_local*lineofsight_ecef;

	/* Calculation of azimuth and elevation angles */
	mp_IrnssSat[sat_id].elevation = -(asin(lineofsight_local(2, 0)));
	mp_IrnssSat[sat_id].azimuth = (atan2(lineofsight_local(1, 0), lineofsight_local(0, 0)));
}