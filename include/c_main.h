#pragma once
#include<string>
#include<random>
#include<Eigen/Dense>
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<math.h>
#include<stdlib.h>
#include<istream>
#include <limits>
#include <vector>

using namespace std;
using namespace Eigen;

class CMain
{
public:
	CMain();
	~CMain();
	int MainFunc();
	void MemoryAlloc(int no_of_iters);
	void MemoryDealloc();
	double *mp_TransTime,
		*mp_IrnssTransTime;

	struct UserMot{

		double  user_pos[3],
				user_vel[3],
				user_acc[3],
				user_jerk[3],
				user_ang[3],
				lat,
				longi,
				height,
				roll,
				pitch,
				yaw;
	};
	UserMot m_UserMotion;

	struct InsOutput {

		Eigen::Matrix<double, 3, 1, Eigen::DontAlign> sf_actual, ang_actual,
		acc_bias, gyro_bias,
		acc_scale, gyro_scale_error,
		acc_misal_error, gyro_misal_error,
		acc_inrun, gyro_inrun,
		acc_noise, gyro_noise;

	};
	InsOutput *m_InsOutput = new InsOutput;
	struct TimeVar{

		unsigned int  week_no,
					days_in_year,
					year;
		double time_of_week;
		double time_of_day;
		int times_rollover;

	};
	TimeVar m_IrnssTimeData;
	TimeVar m_GpsTimeData;

	struct SatData{
		double  x_cord,
				y_cord,
				z_cord,
				x_vel,
				y_vel,
				z_vel;
		double  range,
				range_rt;
		double  azimuth,
			elevation,
			elevation_check;
		double clock_correction;
	};
	
	SatData *mp_GpsSatData;
	SatData *mp_IrnssSatData;

	struct AtmSatDelay{

		double tropo_delay_hop,
		iono_delay_klob;
				//tropo_delay_rtca;
				
	};
	AtmSatDelay m_DelayCalc;
	
	// From Input file
	struct DelayCalcParam{
		double alpha[4],
		beta[4],
		mask_angle,
		humidity;

	};
	DelayCalcParam m_DelayPram;

	void CalcClockError();   //To calculate reciever clock error
	//Random Walk Model
	struct psd{
		long double psd1;
		long double psd2;
		long double dt;
	};
	psd *mp_psd;

	Matrix<long double, 2, 2> m_A;
	Matrix<long double, 2, 2> m_B;

	//Clock error covariance matrix
	Matrix<long double, 2, 2> m_Qclk;
	Matrix<long double, Dynamic, Dynamic> m_Qchol;

	default_random_engine m_rand;
	normal_distribution<double> m_gaussDist;
	Matrix<long double, 2, 1> m_UserErr;

	//To calculate code tracking error
	struct discriminator{
		string type;
		double bandwidth;
		double d;
		double tau_a;
	};
	discriminator m_CodeDisc;
	discriminator m_CarrDisc;
	double elevation_deg;
	double c_n0;

	double f_co; //Code chipping rate
	double f_ca; // Carrier Frequency

	void CalcTrackErr();
	double m_prTrackErr;
	double m_prrTrackErr;

	private:
		void IrnssTimeOfWeek(string time_of_int);
		void GpsTimeOfWeek(string time_of_int);
		void IrnssWeekCrossover();
		void GpsWeekCrossover();
		void ReadUserMotionFile(ifstream& file);
		double  m_TimeOfIter,
			m_TimeStep,
			m_TimeIntoRun,
			m_TimeIntoRunHr,
			m_pseudo_range[31],
			m_pseudo_range_irnss[8];
		static bool end_usermotion_file_flag;

		

};

