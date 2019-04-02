#pragma once
#include<string>
#include "GeographicLib/Geoid.hpp"
using namespace std;
#include "c_main.h"


class CDelayCalc : public CMain
{
public:
	CDelayCalc();
	~CDelayCalc();
	void ComputeDelays(double time_of_week,
						double day_of_year, 
						double lati, 
						double longi,
						double height,
						double azimuth,
						double elevation, 
						struct AtmSatDelay &m_DelayCalc,
						struct DelayCalcParam &DelayPram);

	struct AtmSatDelay m_AtmDelayCalc;

private:
	void IonoDelayKlob();
	void TropoDelayHop();
	void StdAtm(double height_above_sea);


	double m_Lat,
		m_Longi,
		m_Ele,
		m_Azi,
		m_Height,
		m_TimeOfWeek,
		m_DayOfYear;
	
	double tropo_rtca_var;
	double m_Alpha[4],
		m_Beta[4],
		m_Humidity;

	struct MetorologicalParam
	{
		/*
		* Following 5 member variables are meteorological parameters used in RTCA98 model
		* at receiver latitude and day of the year
		*/
		double pressure,
				temp,
				watervap_pr,
				temp_laps_rt,
				watervap_laps_rt;
	};
	MetorologicalParam m_MeteorParam;

	// Pressure, temperature and density for standard atmosphere calculation
	double m_Temp,
		m_Pres,
		m_Density;

	GeographicLib::Geoid geo;
};

