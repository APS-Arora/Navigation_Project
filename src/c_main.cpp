/********************************************************************************************************
* Description            : cpp file of the c_main class, all the member functions are defined here
* Specific library calls : None
* Classes                : 1. CMain
*						   2. CGps
*						   3. CIrnss
*						   4. CDelayCalc
* Assumptions            : None
* Reference              : 1. GPS Interface Control Document IS-GPS-200H
*                          2. Chapter 7, Paul D Grooves : Principles of GNSS, Inertial, and Multisensor Integrated Navigation
*                             System
* Version History        :
* <1.1><Thakur Shivam Singh><20/3/2017>
***********************************************************************************************************/
#include "c_main.h"
#include "c_Gps.h"
#include "c_Irnss.h"
#include "INS.h"
#include "ConstVar.h"
#include "CDelayCalc.h"
#include "COtherSensors.h"

CMain::CMain()
{
	//Allocate memory to structure pointer
	mp_psd = new psd[1];
	//memset(mp_psd, 0, sizeof(psd));

	//Initialize the power spectral density functions
	mp_psd->psd1 = 8.0 * pow(Pi_Const, 2.0) * 2 * pow(10.0, -20.0);
	mp_psd->psd2 = 2.0 * 2.0 * pow(10.0, -19.0);
	mp_psd->dt = (1.0 / 5.0) / pow(10.0, 2.0);

	//Compute the Q matrix
	m_Qclk(0, 0) = mp_psd->psd1*mp_psd->dt + mp_psd->psd2*pow(mp_psd->dt, 3.0) / 3.0;
	m_Qclk(0, 1) = mp_psd->psd2*pow(mp_psd->dt, 2.0) / 2.0;
	m_Qclk(1, 0) = mp_psd->psd2*pow(mp_psd->dt, 2.0) / 2.0;
	m_Qclk(1, 1) = mp_psd->psd2*mp_psd->dt;

	//Cholesky decompostion of Q
	LLT<Matrix<long double, Dynamic, Dynamic>> chol(m_Qclk);
	m_Qchol = chol.matrixL();

	//Random walk model
	m_A << 1.0, mp_psd->dt,
		0.0, 1.0;
	m_B << 1.0, 0.0,
		0.0, 1.0;

	//Initialize the error value 
	m_UserErr(0, 0) = 5.0 * pow(10.0, -5.0);
	m_UserErr(1, 0) = 0.5 * pow(10.0, -8.0);

	m_gaussDist = normal_distribution<double>(0.0, 1.0);

	//Initializing the Code Discriminator 
	m_CodeDisc.bandwidth = 0.5;
	m_CodeDisc.d = 1;
	m_CodeDisc.tau_a = 0.01;
	f_co = 1.023*pow(10, 6);
	m_CodeDisc.type = "DPP"; // Permissable Code Discriminator types: DPP, ELP, ELE, Coh

	//Initializing the Carrier Discriminator 
	m_CarrDisc.bandwidth = 2;
	m_CarrDisc.d = NULL;
	m_CarrDisc.tau_a = 0.01;
	f_ca = 1.023*pow(10, 6);
	m_CarrDisc.type = "FLL"; // Permissable Carrier Discriminator types: Costas, PLL, FLL
}


CMain::~CMain()
{
	delete(mp_psd);
	mp_psd = NULL;

}

/********************************************************************************************************
* Function               : MainFunc
* Description            : Main Function to create objects and call other functions for the main calculations
* Function Parameter     : None
* Return value           : None
* Specific library calls : None
* Functions called       : MemoryAlloc
*						   MemoryDelloc
*						   TimeOfWeek
*						   c_Gps_obj.GpsSat
c_Irnss_obj.IrnssSat
* Assumptions            : GPS and IRNSS Rollover took place on 22 Aug 1999
* Reference              : None
* Version History        :
* <1.1><Thakur Shivam Singh><20/3/2017>
***********************************************************************************************************/

bool CMain::end_usermotion_file_flag = true;

int CMain::MainFunc()
{
	string  strn;
	string time_of_int;

	// Output files
	ofstream ofile_gps("output_gps.txt");
	ofstream ofile_irnss("output_irnss.txt");
	ofstream ofile_gps_delay("gps_delay.txt");
	ofstream ofile_irnss_delay("irnss_delay.txt");
	ofstream ofile_clock_error("receiver_error.txt");
	ofstream ofile_track_error("tracking_error.txt");

	// Output code for INS
	ofstream ofile_ins("ins_output.csv", ios::out | ios::trunc);
	ofile_ins << "Specific_Force_X" << ',' << "Specific_Force_Y" << ',' << "Specific_Force_Z" << ',' << "Angular Rate_X" << ',' << "Angular Rate_Y" << ',' << "Angular Rate_Z" << '\n';

	// Output File for Barometric Altimeter and Magnetometer
	ofstream ofile_osensors("sensed_output_other_sensors.csv");
	ofile_osensors << "Time, Magnetic Reading X, Magnetic Reading Y, Magnetic Reading Z , Magnetic Heading, Pressure\n";
	ofstream ofile_osensors_truth("true_output_other_sensors.csv");
	ofile_osensors_truth << "Time, Magnetic Reading X, Magnetic Reading Y, Magnetic Reading Z , Magnetic Heading, Pressure\n";

	// Precision for output files
	ofile_gps.precision(15);
	ofile_irnss.precision(15);
	ofile_gps_delay.precision(15);
	ofile_irnss_delay.precision(15);
	ofile_ins.precision(15);
	ofile_osensors.precision(15);
	ofile_osensors_truth.precision(15);
	std::cout.precision(10);
	/* Input file */
	ifstream iFile("UserConfigInput.txt");
	ifstream  iFile_1("UserMotionInput.csv");

	// Ignore first 2 lines of User Motion Input File
	iFile_1.ignore(numeric_limits<streamsize>::max(), '\n');
	iFile_1.ignore(numeric_limits<streamsize>::max(), '\n');
	/* Reading input from the UserConfigInput.txt file */
	string strin;
	double doubl;
	int inte;
	while (!iFile.eof())
	{

		iFile >> strin; iFile >> strin; iFile >> strin; iFile >> strin;
		time_of_int = strin;
		iFile >> strin; iFile >> strin; iFile >> doubl;
		m_TimeStep = doubl;
		iFile >> strin; iFile >> strin; iFile >> strin; iFile >> doubl;
		m_TimeOfIter = doubl;
		iFile >> strin; iFile >> strin; iFile >> doubl;
		m_DelayPram.mask_angle = doubl;
		iFile >> strin; iFile >> strin; iFile >> doubl; 
		m_DelayPram.alpha[0] = doubl;
		iFile >> doubl;
		m_DelayPram.alpha[1] = doubl;
		iFile >> doubl;
		m_DelayPram.alpha[2] = doubl;
		iFile >> doubl;
		m_DelayPram.alpha[3] = doubl;
		iFile >> strin; iFile >> strin; iFile >> doubl;
		m_DelayPram.beta[0] = doubl;
		iFile >> doubl;
		m_DelayPram.beta[1] = doubl;
		iFile >> doubl;
		m_DelayPram.beta[2] = doubl;
		iFile >> doubl;
		m_DelayPram.beta[3] = doubl;
		iFile >> strin; iFile >> doubl;
		m_DelayPram.humidity = doubl;
		/*iFile >> strin; iFile >> strin; iFile >> strin; iFile >> inte;
		m_DelayPram.tropo_flag = inte;*/

	}

	// Calculate the number of times simulation will be done in a given time of simulation
	int no_of_iters;
	no_of_iters = int(m_TimeOfIter / m_TimeStep) + 1;

	// Function call to allocate memory to the struct pointers of main class
	MemoryAlloc(no_of_iters);

	// Function call to convert time in date-time to GPS 'time of week (sec)' and 'week number'
	IrnssTimeOfWeek(time_of_int);
	GpsTimeOfWeek(time_of_int);

	// Create object for delay estimation
	CDelayCalc c_Delay_obj;

	// Create GPS object
	CGps c_Gps_obj;

	// Create IRNSS object
	CIrnss c_Irnss_obj;

	//Create INS object
	INS ins_obj;

	// Create object for Other Sensors
	COtherSensors os_obj("OtherSensorsConfig.csv");

	int time_iter = 0;

	// For loop to calculate all variables for each time step
	for (double itervar = 0; itervar <= m_TimeOfIter; itervar += m_TimeStep)
	{

		// Iteration number
		int iter_no = int(itervar / m_TimeStep);

		m_TimeIntoRun = itervar;
		//int csv_row_no = iter_no + 4;
		int csv_row_no = time_iter + 4;

		// Function call to read user motion file
		if (end_usermotion_file_flag == true)
		{
			ReadUserMotionFile(iFile_1);
		}


		//Increase the time of week by the time step for each iteration
		//m_IrnssTimeData.time_of_week = m_IrnssTimeData.time_of_week + itervar;

		
		//Function to account for the week crossover for IRNSS
		if (m_IrnssTimeData.time_of_week > 604800)
		{
			m_IrnssTimeData.week_no = m_IrnssTimeData.week_no + 1;
			m_IrnssTimeData.time_of_week = m_IrnssTimeData.time_of_week - 604800;

			if (m_IrnssTimeData.week_no > 1023)
			{
				m_IrnssTimeData.week_no = m_IrnssTimeData.week_no - 1024;
				m_IrnssTimeData.times_rollover = m_IrnssTimeData.times_rollover + 1;
			}
		}
		
		
		//Function to account for the week crossover for GPS
		if (m_GpsTimeData.time_of_week > 604800)
		{
			m_GpsTimeData.week_no = m_GpsTimeData.week_no + 1;
			m_GpsTimeData.time_of_week = m_GpsTimeData.time_of_week - 604800;

			if (m_GpsTimeData.week_no > 1023)
			{
				m_GpsTimeData.week_no = m_GpsTimeData.week_no - 1024;
				m_GpsTimeData.times_rollover = m_GpsTimeData.times_rollover + 1;
			}
		}

		// Change time in hour so as to save in the output file
		m_TimeIntoRunHr = itervar / 3600.0;

		//Calling INS main function
		ins_obj.main_function(m_UserMotion, m_InsOutput);

		char comma = ',';

		//Writing the outputs to file
		ofile_ins << m_InsOutput->sf_actual[0] << comma << m_InsOutput->sf_actual[1] << comma << m_InsOutput->sf_actual[2] << comma << m_InsOutput->ang_actual[0] << comma << m_InsOutput->ang_actual[1] << comma << m_InsOutput->ang_actual[2] << comma;
		ofile_ins << m_InsOutput->acc_bias[0] << comma << m_InsOutput->acc_bias[1] << comma << m_InsOutput->acc_bias[2] << comma << m_InsOutput->gyro_bias[0] << comma << m_InsOutput->gyro_bias[1] << comma << m_InsOutput->gyro_bias[2] << comma;
		ofile_ins << m_InsOutput->acc_scale[0] << comma << m_InsOutput->acc_scale[1] << comma << m_InsOutput->acc_scale[2] << comma << m_InsOutput->gyro_scale_error[0] << comma << m_InsOutput->gyro_scale_error[1] << comma << m_InsOutput->gyro_scale_error[2] << comma;
		ofile_ins << m_InsOutput->acc_misal_error[0] << comma << m_InsOutput->acc_misal_error[1] << comma << m_InsOutput->acc_misal_error[2] << comma << m_InsOutput->gyro_misal_error[0] << comma << m_InsOutput->gyro_misal_error[1] << comma << m_InsOutput->gyro_misal_error[2] << comma;
		ofile_ins << m_InsOutput->acc_inrun[0] << comma << m_InsOutput->acc_inrun[1] << comma << m_InsOutput->acc_inrun[2] << comma << m_InsOutput->gyro_inrun[0] << comma << m_InsOutput->gyro_inrun[1] << comma << m_InsOutput->gyro_inrun[2] << comma;
		ofile_ins << m_InsOutput->acc_noise[0] << comma << m_InsOutput->acc_noise[1] << comma << m_InsOutput->acc_noise[2] << comma << m_InsOutput->gyro_noise[0] << comma << m_InsOutput->gyro_noise[1] << comma << m_InsOutput->gyro_noise[2] << '\n';

		// Taking Altimeter and Magnetic Measurements
		os_obj.m_Measure(m_GpsTimeData, m_UserMotion);

		//Writing the outputs to file
		ofile_osensors << itervar << ',' << os_obj.m_RecentMagMeasurement.mag_strength.x() << ',' << os_obj.m_RecentMagMeasurement.mag_strength.y() << ',' << os_obj.m_RecentMagMeasurement.mag_strength.z() << ',' << os_obj.m_RecentMagMeasurement.mag_heading << ',' << os_obj.m_RecentPressMeasurement << '\n';
		ofile_osensors_truth << itervar << ',' << os_obj.m_RecentMagTruth.mag_strength.x() << ',' << os_obj.m_RecentMagTruth.mag_strength.y() << ',' << os_obj.m_RecentMagTruth.mag_strength.z() << ',' << os_obj.m_RecentMagTruth.mag_heading << ',' << os_obj.m_RecentPressTruth << '\n';

		//Calculate reciever clock error
		CalcClockError();

		//Write output of calculated clock error
		ofile_clock_error << m_UserErr(0, 0) << " ";
		ofile_clock_error << m_UserErr(1, 0) << endl;

		// Function call from GPS class
		c_Gps_obj.GpsSat(m_GpsTimeData, m_UserMotion, mp_GpsSatData, mp_TransTime);
		
		// Pass the value of calculated GPS satellite variables back to the main class
		for (int gps_sat_id = 0; gps_sat_id < MaxGpsSat_Const; gps_sat_id++)
		{

			c_Delay_obj.ComputeDelays(m_GpsTimeData.time_of_week,
				m_GpsTimeData.days_in_year,
				m_UserMotion.lat,
				m_UserMotion.longi,
				m_UserMotion.height,
				mp_GpsSatData[gps_sat_id].azimuth,
				mp_GpsSatData[gps_sat_id].elevation,
				m_DelayCalc,
				m_DelayPram);

			//Calculate code tracking error
			CalcTrackErr(abs(mp_GpsSatData[gps_sat_id].elevation)*180.0 / Pi_Const);

			// Calculate pseudo range and pseudo range rate		
				m_pseudo_range[gps_sat_id] = mp_GpsSatData[gps_sat_id].range + m_DelayCalc.iono_delay_klob + m_DelayCalc.tropo_delay_hop -
					mp_GpsSatData[gps_sat_id].clock_correction + SpeedLight_Const*m_UserErr(0, 0) + m_prTrackErr;
				m_pseudo_range_rate[gps_sat_id] = mp_GpsSatData[gps_sat_id].range_rt + SpeedLight_Const*m_UserErr(1, 0) + m_prrTrackErr;

			// Write output for the GPS data
			ofile_gps << mp_GpsSatData[gps_sat_id].x_cord << "  ";
			ofile_gps << mp_GpsSatData[gps_sat_id].y_cord << "  ";
			ofile_gps << mp_GpsSatData[gps_sat_id].z_cord << "  ";
			ofile_gps << mp_GpsSatData[gps_sat_id].x_vel << "  ";
			ofile_gps << mp_GpsSatData[gps_sat_id].y_vel << "  ";
			ofile_gps << mp_GpsSatData[gps_sat_id].z_vel << "  ";
			ofile_gps << mp_GpsSatData[gps_sat_id].range << "  ";
			ofile_gps << mp_GpsSatData[gps_sat_id].range_rt << "  ";
			ofile_gps << mp_GpsSatData[gps_sat_id].azimuth << "  ";
			ofile_gps << mp_GpsSatData[gps_sat_id].elevation << "  ";
			ofile_gps << m_TimeIntoRunHr << "  ";
			ofile_gps << m_pseudo_range[gps_sat_id] << "  ";
			ofile_gps << m_pseudo_range_rate[gps_sat_id] << "  ";
			ofile_gps << mp_GpsSatData[gps_sat_id].clock_correction << endl;


			//Write outbput for calculated delay
			ofile_gps_delay << m_DelayCalc.iono_delay_klob << " ";
			ofile_gps_delay << m_DelayCalc.tropo_delay_hop << " ";
			ofile_gps_delay << mp_GpsSatData[gps_sat_id].azimuth << " ";
			ofile_gps_delay << mp_GpsSatData[gps_sat_id].elevation << "  ";
			ofile_gps_delay << m_TimeIntoRunHr << endl;
			
			//Write tracking error to file
			ofile_track_error << m_prTrackErr << " " << m_prrTrackErr << " ";
		}
	
		
		//Function call from IRNSS class
		c_Irnss_obj.IrnssSat(m_IrnssTimeData, m_UserMotion, mp_IrnssSatData, mp_IrnssTransTime);

		// Write output for the IRNSS data
		for (int irnss_sat_id = 0; irnss_sat_id < MaxIrnssSat_Const; irnss_sat_id++)
		{

			c_Delay_obj.ComputeDelays(m_IrnssTimeData.time_of_week,
										m_IrnssTimeData.days_in_year,
										m_UserMotion.lat,
										m_UserMotion.longi,
										m_UserMotion.height,
										mp_IrnssSatData[irnss_sat_id].azimuth,
										mp_IrnssSatData[irnss_sat_id].elevation,
										m_DelayCalc,
										m_DelayPram);

			//Calculate code tracking error
			CalcTrackErr(abs(mp_IrnssSatData[irnss_sat_id].elevation)*180.0 / Pi_Const);

			// Calculate pseudo range after adding iono and tropo delay and subtracting clock error			
			m_pseudo_range_irnss[irnss_sat_id] = mp_IrnssSatData[irnss_sat_id].range + m_DelayCalc.iono_delay_klob + m_DelayCalc.tropo_delay_hop -
				mp_IrnssSatData[irnss_sat_id].clock_correction + SpeedLight_Const*m_UserErr(0, 0) + m_prTrackErr;
			m_pseudo_range_rate_irnss[irnss_sat_id] = mp_IrnssSatData[irnss_sat_id].range_rt + SpeedLight_Const*m_UserErr(1, 0) + m_prrTrackErr;


			ofile_irnss << mp_IrnssSatData[irnss_sat_id].x_cord << "  ";
			ofile_irnss << mp_IrnssSatData[irnss_sat_id].y_cord << "  ";
			ofile_irnss << mp_IrnssSatData[irnss_sat_id].z_cord << "  ";
			ofile_irnss << mp_IrnssSatData[irnss_sat_id].x_vel << "  ";
			ofile_irnss << mp_IrnssSatData[irnss_sat_id].y_vel << "  ";
			ofile_irnss << mp_IrnssSatData[irnss_sat_id].z_vel << "  ";
			ofile_irnss << mp_IrnssSatData[irnss_sat_id].range << "  ";
			ofile_irnss << mp_IrnssSatData[irnss_sat_id].range_rt << "  ";
			ofile_irnss << mp_IrnssSatData[irnss_sat_id].azimuth << "  ";
			ofile_irnss << mp_IrnssSatData[irnss_sat_id].elevation << "  ";
			ofile_irnss << m_TimeIntoRunHr << "  ";
			ofile_irnss << m_pseudo_range_irnss[irnss_sat_id] << "  ";
			ofile_irnss << m_pseudo_range_rate_irnss[irnss_sat_id] << "  ";
			ofile_irnss << mp_IrnssSatData[irnss_sat_id].clock_correction << endl;

			//Write outbput for calculated delay
			ofile_irnss_delay << m_DelayCalc.iono_delay_klob << " ";
			ofile_irnss_delay << m_DelayCalc.tropo_delay_hop << " ";
			ofile_irnss_delay << mp_IrnssSatData[irnss_sat_id].azimuth << " ";
			ofile_irnss_delay << mp_IrnssSatData[irnss_sat_id].elevation << "  ";
			ofile_irnss_delay << m_TimeIntoRunHr << endl;

			//Write tracking error to file
			ofile_track_error << m_prTrackErr << " " << m_prrTrackErr << endl;
		}
		
		

		//Increase the time of week by the time step after each iteration
		m_GpsTimeData.time_of_week = m_GpsTimeData.time_of_week + m_TimeStep;

		// Increase the time of current day by the time step after each iteration
		m_GpsTimeData.time_of_day = m_GpsTimeData.time_of_day + m_TimeStep;
		if (m_GpsTimeData.time_of_day >= 86400)
		{
			m_GpsTimeData.time_of_day = m_GpsTimeData.time_of_day - 86400;
			m_GpsTimeData.days_in_year += 1;
			if (m_GpsTimeData.days_in_year >= 366 && (m_GpsTimeData.year % 4) == 0)
			{
				m_GpsTimeData.days_in_year -= 366;
				m_GpsTimeData.year += 1;
			}
			if (m_GpsTimeData.days_in_year >= 365 && (m_GpsTimeData.year % 4) != 0)
			{
				m_GpsTimeData.days_in_year -= 365;
				m_GpsTimeData.year += 1;
			}

		}

		
		// For irnss satellites
		//Increase the time of week by the time step after each iteration
		m_IrnssTimeData.time_of_week = m_IrnssTimeData.time_of_week + m_TimeStep;

		// Increase the time of current day by the time step after each iteration
		m_IrnssTimeData.time_of_day = m_IrnssTimeData.time_of_day + m_TimeStep;
		if (m_IrnssTimeData.time_of_day > 86400)
		{
			m_IrnssTimeData.time_of_day = m_IrnssTimeData.time_of_day - 86400;
			m_IrnssTimeData.days_in_year += 1;
			if (m_IrnssTimeData.days_in_year >= 366 && (m_IrnssTimeData.year % 4) == 0)
			{
				m_IrnssTimeData.days_in_year -= 366;
				m_IrnssTimeData.year += 1;
			}
			if (m_IrnssTimeData.days_in_year >= 365 && (m_IrnssTimeData.year % 4) != 0)
			{
				m_IrnssTimeData.days_in_year -= 365;
				m_IrnssTimeData.year += 1;
			}
		}
		


		std::cout << "iter_no " << iter_no << endl;
		//cout << "time_iter " << time_iter << endl;
		//cout << "time_of_week " << m_GpsTimeData.time_of_week << endl;
		time_iter++;
		
	}


	// Function call to deallocate the memory of struct pointers of main class
	MemoryDealloc();

	std::cout << "no_of_iters" << no_of_iters << endl;
	
	std::system("pause");
	return 0;
}


/********************************************************************************************************
* Function               : IrnssTimeOfWeek
* Description            : Function to calculate time of week and week number for the irnss consellation from the date-time format
* Function Parameter     : time of interest: Input time in date time format
* Return value           : None
* Specific library calls : None
* Functions called       : None
* Assumptions            : GPS and IRNSS Rollover took place on 22 Aug 1999
* Reference              : None
* Version History        :
* <1.1><Thakur Shivam Singh><20/3/2017>
***********************************************************************************************************/

void CMain::IrnssTimeOfWeek(string time_of_int){

	int year, month, day, hour, min, sec;
	int days_in_year, days_aft_rlovr;

	/* Use substrings of string time_of_int */
	year = stoi(time_of_int.substr(0, 4));
	month = stoi(time_of_int.substr(5, 7));
	day = stoi(time_of_int.substr(8, 10));
	hour = stoi(time_of_int.substr(11, 13));
	min = stoi(time_of_int.substr(14, 16));
	sec = stoi(time_of_int.substr(17, 19));

	m_IrnssTimeData.times_rollover = 0;
	m_IrnssTimeData.year = year;
	// Calculate number of seconds passed in the current day
	m_IrnssTimeData.time_of_day = (3600 * hour) + (60 * min) + sec;


	/* Calculate no of days passed in the year till the given date  */
	/* If current year is a leap year */
	if (year % 4 == 0)
	{

		if ((month - 2) / 6 != 0)
		{
			days_in_year = (month - 1) * 30 + ((month - 1) / 2 + 1) - 1 + day;
		}
		else
		{
			if (month < 3)
			{
				days_in_year = (month - 1) * 31 + day;
			}
			else
			{
				days_in_year = (month - 1) * 30 + ((month - 1) / 2) - 1 + day;
			}
		}
	}

	/* If current year is not a leap year */
	else
	{
		if ((month - 2) / 6 != 0){
			days_in_year = (month - 1) * 30 + ((month - 1) / 2 + 1) - 2 + day;
		}
		else
		{

			if (month < 3)
			{
				days_in_year = (month - 1) * 31 + day;
			}
			else
			{
				days_in_year = (month - 1) * 30 + ((month - 1) / 2) - 2 + day;
			}
		}
	}

	m_IrnssTimeData.days_in_year = days_in_year;

	/* Calulate no of days between time of IRNSS epoch(22 Aug 1999) and Time of interest */
	if (year > 2000)
	{
		// days_aft_rlovr = 131 + ((stoi(time_of_int.substr(2, 2))) * 365)
		//	+ ((stoi(time_of_int.substr(2, 2))) / 4 + 1) + days_in_year;

		days_aft_rlovr = 131 + ((year - 2000) * 365)
				+ (((year - 2000) / 4) + 1) + days_in_year;
	}

	if (year == 2000)
	{
		days_aft_rlovr = 131 + days_in_year;
	}

	if (year < 2000)
	{
		days_aft_rlovr = days_in_year - 131;
	}

	/* Calulate no of weeks between time of GPS & IRNSS rollover(22 Aug 1999) and Time of interest */
	m_IrnssTimeData.week_no = days_aft_rlovr / 7;

	while (m_IrnssTimeData.week_no > 1023)
	{
		m_IrnssTimeData.week_no = m_IrnssTimeData.week_no - 1023;
		m_IrnssTimeData.times_rollover = m_IrnssTimeData.times_rollover + 1;
		days_aft_rlovr = days_aft_rlovr - (1024 * 7);
	}

	/* No of seconds passed in the current week */
	m_IrnssTimeData.time_of_week = (days_aft_rlovr - m_IrnssTimeData.week_no * 7) * 24 * 60 * 60
							+ hour * 60 * 60 + min * 60 + sec;
}

/********************************************************************************************************
* Function               : GpsTimeOfWeek
* Description            : Function to deallocate memory to struct pointers in main class
* Function Parameter     : None
* Return value           : None
* Specific library calls : None
* Functions called       : None
* Assumptions            : None
* Reference              : None
* Version History        :
* <1.1><Thakur Shivam Singh><20/3/2017>
***********************************************************************************************************/

void CMain::GpsTimeOfWeek(string time_of_int){

	int year, month, day, hour, min, sec;
	int days_in_year, days_aft_rlovr;

	/* Use substrings of string time_of_int */
	year = stoi(time_of_int.substr(0, 4));
	month = stoi(time_of_int.substr(5, 7));
	day = stoi(time_of_int.substr(8, 10));
	hour = stoi(time_of_int.substr(11, 13));
	min = stoi(time_of_int.substr(14, 16));
	sec = stoi(time_of_int.substr(17, 19));

	m_GpsTimeData.year = year;
	m_GpsTimeData.time_of_day = (3600 * hour) + (60 * min) + sec;

	/* Calculate no of days passed in the year till the given date  */
	/* If current year is a leap year */
	if (year % 4 == 0)
	{

		if ((month - 2) / 6 != 0)
		{
			days_in_year = (month - 1) * 30 + ((month - 1) / 2 + 1) - 1 + day;
		}
		else
		{
			if (month < 3)
			{
				days_in_year = (month - 1) * 31 + day;
			}
			else
			{
				days_in_year = (month - 1) * 30 + ((month - 1) / 2) - 1 + day;
			}
		}
	}

	/* If current year is not a leap year */
	else
	{
		if ((month - 2) / 6 != 0){
			days_in_year = (month - 1) * 30 + ((month - 1) / 2 + 1) - 2 + day;
		}
		else
		{

			if (month < 3)
			{
				days_in_year = (month - 1) * 31 + day;
			}
			else
			{
				days_in_year = (month - 1) * 30 + ((month - 1) / 2) - 2 + day;
			}
		}
	}

	m_GpsTimeData.days_in_year = days_in_year;

	/* Calulate no of days between time of GPS epoch (6 Jan 1980) and Time of interest */
	if (year > 1980)
	{
		// days_aft_rlovr = 131 + ((stoi(time_of_int.substr(2, 2))) * 365)
		//	+ ((stoi(time_of_int.substr(2, 2))) / 4 + 1) + days_in_year;

		days_aft_rlovr = 360 + ((year - 1981) * 365)
			+ (((year - 1981) / 4)) + days_in_year;
	}

	if (year == 1980)
	{
		days_aft_rlovr = days_in_year - 6;
	}

	/* Calulate no of weeks between time of GPS & IRNSS rollover(22 Aug 1999) and Time of interest */
	m_GpsTimeData.week_no = days_aft_rlovr / 7;

	while (m_GpsTimeData.week_no > 1023)
	{
		m_GpsTimeData.week_no = m_GpsTimeData.week_no - 1024;
		m_GpsTimeData.times_rollover = m_GpsTimeData.times_rollover + 1;
		days_aft_rlovr = days_aft_rlovr - (1024 * 7);
	}

	/* No of seconds passed in the current week */
	m_GpsTimeData.time_of_week = (days_aft_rlovr - m_GpsTimeData.week_no * 7) * 24 * 60 * 60
		+ hour * 60 * 60 + min * 60 + sec;
}

/********************************************************************************************************
* Function               : ReadUserMotionFile
* Description            : Function to read user motion file
* Function Parameter     : file: Input file
*						   row: row number of input file
* Return value           : None
* Specific library calls : None
* Functions called       : None
* Assumptions            : None
* Reference              : None
* Version History        :
* <1.1><Thakur Shivam Singh><20/3/2017>
***********************************************************************************************************/
void CMain::ReadUserMotionFile(ifstream& file)
{
	string param;
	double time;

	if (!getline(file, param, ','))
	{
		m_UserMotion.user_vel[0] = 0;
		m_UserMotion.user_vel[1] = 0;
		m_UserMotion.user_vel[2] = 0;

		m_UserMotion.user_acc[0] = 0;
		m_UserMotion.user_acc[1] = 0;
		m_UserMotion.user_acc[2] = 0;

		m_UserMotion.user_jerk[0] = 0;
		m_UserMotion.user_jerk[1] = 0;
		m_UserMotion.user_jerk[2] = 0;

		end_usermotion_file_flag = false;
	}

	if (getline(file, param, ','))
	{
		//Reading time into run
		//getline(file, param, ',');
		//time = stod(param);

		// Reading user position
		m_UserMotion.user_pos[0] = stod(param);
		getline(file, param, ',');
		m_UserMotion.user_pos[1] = stod(param);
		getline(file, param, ',');
		m_UserMotion.user_pos[2] = stod(param);

		// Read user velocity
		getline(file, param, ',');
		m_UserMotion.user_vel[0] = stod(param);
		getline(file, param, ',');
		m_UserMotion.user_vel[1] = stod(param);
		getline(file, param, ',');
		m_UserMotion.user_vel[2] = stod(param);

		// Read user acceleration
		getline(file, param, ',');
		m_UserMotion.user_acc[0];
		getline(file, param, ',');
		m_UserMotion.user_acc[1] = stod(param);
		getline(file, param, ',');
		m_UserMotion.user_acc[2] = stod(param);

		// Read user Jerk
		getline(file, param, ',');
		m_UserMotion.user_jerk[0] = stod(param);
		getline(file, param, ',');
		m_UserMotion.user_jerk[1] = stod(param);
		getline(file, param, ',');
		m_UserMotion.user_jerk[2] = stod(param);

		//Read user Latitude, Longitude and height
		getline(file, param, ',');
		m_UserMotion.lat = stod(param);
		getline(file, param, ',');
		m_UserMotion.longi = stod(param);
		getline(file, param, ',');
		m_UserMotion.height = stod(param);

		//Read user Orientation
		getline(file, param, ',');
		m_UserMotion.yaw = stod(param);
		getline(file, param, ',');
		m_UserMotion.pitch = stod(param);
		getline(file, param, ',');
		m_UserMotion.roll = stod(param);

		//Read user Angular Velocity
		getline(file, param, ',');
		m_UserMotion.user_ang[0] = stod(param);
		getline(file, param, ',');
		m_UserMotion.user_ang[1] = stod(param);
		getline(file, param, ',');
		m_UserMotion.user_ang[2] = stod(param);

		getline(file, param);
	}

}

/********************************************************************************************************
* Function               : MemoryAlloc
* Description            : Function to allocate memory to struct pointers in main class
* Function Parameter     : None
* Return value           : None
* Specific library calls : None
* Functions called       : None
* Assumptions            : None
* Reference              : None
* Version History        :
* <1.1><Thakur Shivam Singh><20/3/2017>
***********************************************************************************************************/

void CMain::MemoryAlloc(int no_of_iters)
{

	mp_GpsSatData = new SatData[MaxGpsSat_Const];
	memset(mp_GpsSatData, 0, sizeof(SatData)*MaxGpsSat_Const);

	mp_TransTime = new double[MaxGpsSat_Const];
	memset(mp_TransTime, 0, sizeof(double)*MaxGpsSat_Const);


	mp_IrnssSatData = new SatData[MaxIrnssSat_Const];
	memset(mp_IrnssSatData, 0, sizeof(SatData)*MaxIrnssSat_Const);

	mp_IrnssTransTime = new double[MaxIrnssSat_Const];
	memset(mp_IrnssTransTime, 0, sizeof(double)*MaxIrnssSat_Const);

	/*
	mp_gps_sat_data = new SatData*[MaxGpsSat_Const];

	for (int iter = 0; iter < MaxGpsSat_Const; iter++)
	{
		mp_gps_sat_data[iter] = new SatData[no_of_iters];
		memset(mp_gps_sat_data[iter], 0, no_of_iters*sizeof(SatData));
	}


	
	mp_irnss_sat_data = new SatData*[MaxIrnssSat_Const];
	for (int jter = 0; jter < MaxIrnssSat_Const; jter++)
	{
		mp_irnss_sat_data[jter] = new SatData[no_of_iters];
		memset(mp_irnss_sat_data[jter], 0, sizeof(SatData)*no_of_iters);
	}
	*/

}


/********************************************************************************************************
* Function               : MemoryDealloc
* Description            : Function to deallocate memory to struct pointers in main class
* Function Parameter     : None
* Return value           : None
* Specific library calls : None
* Functions called       : None
* Assumptions            : None
* Reference              : None
* Version History        :
* <1.1><Thakur Shivam Singh><20/3/2017>
***********************************************************************************************************/

void CMain::MemoryDealloc()
{

	delete(mp_GpsSatData);
	mp_GpsSatData = NULL;

	delete(mp_TransTime);
	mp_TransTime = NULL;

	delete(mp_IrnssSatData);
	mp_IrnssSatData = NULL;
	/*
	for (int iter = 0; iter < MaxGpsSat_Const; iter++)
	{
		delete[] mp_gps_sat_data[iter];
	}


	delete[] mp_gps_sat_data;
	mp_gps_sat_data = NULL;
	
	for (int iter = 0; iter < MaxIrnssSat_Const; iter++)
	{
		delete[] mp_irnss_sat_data[iter];
	}

	delete[] mp_irnss_sat_data;
	mp_irnss_sat_data = NULL;
	*/
}


/********************************************************************************************************
* Function               : CalcClockError
* Description            : Function to calculate reciever clock bias and drift
* Function Parameter     : None
* Return value           : None
* Specific library calls : #include<random>
* Functions called       : None
* Assumptions            : None
* Reference              : Random Walk model
* Version History        :
* 
***********************************************************************************************************/

void CMain::CalcClockError(){
	long int n_iter = 0.1 / mp_psd->dt;

	Matrix<long double, 2, 1> wn_bias_drift;
	Matrix<long double, 2, 1> wnc;


	for (long int i = 0; i < n_iter; i++){
		wn_bias_drift << m_gaussDist(m_rand),
			m_gaussDist(m_rand);
		wnc = m_Qchol*wn_bias_drift;

		m_UserErr = m_A*m_UserErr + m_B*wnc;

	}
}

/********************************************************************************************************
* Function               : CalcTrackErr
* Description            : Function to calculate code tracking error
* Function Parameter     : None
* Return value           : None
* Specific library calls : #include<random>
* Functions called       : None
* Assumptions            : None
* Reference              : Ch.7 Paul D Groves
* Version History        :
*
***********************************************************************************************************/

void CMain::CalcTrackErr(double elevation_deg){
	
	double sdTrackErr;
	c_n0 = pow(10, 2.0) + ((pow(10, 4.5) - pow(10, 2.0)) / (75 - 5))*(elevation_deg-5);

	if (m_CodeDisc.type == "DPP")
		//Standard Deviation for DPP
		sdTrackErr = (SpeedLight_Const / f_co)*sqrt((m_CodeDisc.bandwidth*m_CodeDisc.d / 2 / c_n0)*(1 + 1 / c_n0 / m_CodeDisc.tau_a));
	else if (m_CodeDisc.type == "ELP")
		//Standard Deviation for ELP
		sdTrackErr = (SpeedLight_Const / f_co)*sqrt((m_CodeDisc.bandwidth*m_CodeDisc.d / 2 / c_n0)*(1 + 2 / (2 - m_CodeDisc.d) / c_n0 / m_CodeDisc.tau_a));
	else if (m_CodeDisc.type == "ELE" || m_CodeDisc.type == "Coh")
		//Standard Deviation for ELE or Coh
		sdTrackErr = (SpeedLight_Const / f_co)*sqrt(m_CodeDisc.bandwidth*m_CodeDisc.d / 2 / c_n0);

	m_prTrackErr = sdTrackErr*m_gaussDist(m_rand);

	if (m_CarrDisc.type == "FLL")
		sdTrackErr = SpeedLight_Const*sqrt(4 * m_CarrDisc.bandwidth*(1 + 1 / (c_n0*m_CarrDisc.tau_a)) / c_n0) / (2 * EIGEN_PI*f_ca*m_CarrDisc.tau_a);

	m_prrTrackErr = sdTrackErr*m_gaussDist(m_rand);
}