
/*
******************************************************************************
*  File                :  GnssModelConsts.h
*
*  Description         :  This header file contains the definition
*                         of all constants used in GNSSModel dll
*
*  Compiler            :
*
*
*  Compiler options    :
*
*
*  H/W platform        :  NA
*
*  Portability         :  None
*
*  Author(s)           :
*
*  Classes             :  CGnssDopGen
*                         CGnssAtmosphericModel
*                         CGnssRecAntennaPattern
*                         CGnssUserMotionModel
*
*  References          :
*
*  Version History     :
*  <Version Number>   <Author>     <date>    <defect Number>  <Modification
*                                                             made and the
*                                                             reason for
*                                                             modification>
*
*********************************************************************************************/

const double GMProd_Const = 3.986004418E+14,
			 Pi_Const = 3.141592653589793,
			 RtEarthRotn_Const = 7.2921151467 / 100000,
			 SecsInWeek_Const = 604800,
			 SpeedLight_Const = 299792458,
			 PolarRadiusOfEarth_Const = 6356752.3142,
			 EqtRadiusOfEarth_Const = 6378137.0,
			 RadiusOfEarth_Const = 6356766.0,
			 Exp_Const = 2.7183,
			 Eccentricity_Const = 0.0818191908425;

const int MaxGpsSat_Const = 31,
		MaxIrnssSat_Const = 7,
		Limit_Ele_Const = 0;


// Integer constants
const int SixtyFiveConst = 65;
const int SixtyEightConst = 68;
const int SeventyThreeConst = 73;
const int SeventyFiveConst = 75;
const int EightyConst = 80;
const int OneNotFiveConst = 105;
const int OneHundredSixConst = 106;
const int OneTenConst = 110;
const int TwoHundredOneConst = 201;
const int ThreeThousandConst = 3000;
const int SixEightyFourConst = 684;

// Number of seconds, minutes, hours, etc.
const double DaysInYear_Const = 365.25;
const int DaysFrmMarToDecConst = 306;

const double NumOfSecIn20HrsConst = 72000;
const double NumOfSecIn14HrsConst = 50400;
const int NumOfMonthsInYrConst = 12;
const int DaysInOneMonthConst = 30;


// Constants for mean sea level height computation
const int MSeaTableGrdSizConst = 10;
const int GeoidTabLatRowsConst = 36;


// All trpospheric model related constants
const double TropoMaxHeight_Const = 25000;
const int Tropo1255Const = 1255;
const double MSeaCorrLowLimConst = -106;
const double MSeaCorrUppLimConst = 76;
const double UniversalGasConst = 8.31447;
const double TwoThreeSevPt3Const = (double)237.3;
const double SevenPtFiveConst = (double)7.5;
const double MolMassDryAirConst = 0.0289644;
const double DegToKelvinConst = (double)273.15;


const int MinDayOfYrNorth_Const = 28;
const int MinDayOfYrSouth_Const = 211;
const int MaxWetHeightConst = 11000;
const double IppGridHeightConst = 350000;
const double EarthHeightConst = 350000;
const double ZeroPtZeroFiveConst = 0.05;
const double PtZeroOneFiveConst = 0.015;
const double ZeroPtOneOneConst = 0.11;
const double PtTwoFiveConst = 0.25;
const double ThreePtSevenFiveConst = 3.75;

const double ZeroPtFiveThreeConst = 0.53;
const double PtZeroTwoTwoConst = 0.022;
const double PtFourOneSixConst = 0.416;
const double OnePtSixOneSevConst = 1.617;
const double PtZeroSixFourConst = 0.064;
const double KlobucharThConst = 1.57;//1.59245043403625;
const double TwoPtTwoFiveConst = 2.25;
const double FiveExpMinus9Const = 0.000000005;
const double SixPtTwoFiveConst = 6.25;
const double SevenNinePtFiveConst = 79.5;
const double OneExpMinusSixConst = 1e-6;

// Constants used for water vapor pressure calculation
const double SatPressureConst = (double)6.1078;

// Constants used for meteorological parameters computation used in all tropospheric models
// Average values of dry pressure
const double TropoPressure_Const[5] = { 1013.25, 1017.25, 1015.75, 1011.75, 1013.00 };

// Average values of temperature
const double TropoTemp_Const[5] = { 299.65, 294.15, 283.15, 272.15, 263.65 };

// Average values of water vapor pressure
const double TropoWaterVapPressure_Const[5] = { 26.31, 21.79, 11.66, 6.78, 4.11 };

// Average values of temperature lapse rate
const double TempLapsRt_Const[5] = { 6.30e-3, 6.05e-3, 5.58e-3, 5.39e-3, 4.53e-3 };

// Average values of water vapor pressure lapse rate
const double WaterVapPressureLapsRt_Const[5] = { 2.77, 3.15, 2.57, 1.81, 1.55 };

// Seasonal variation of dry pressure
const double DeltaPressure_Const[5] = { 0, -3.75, -2.25, -1.75, -0.50 };

// Seasonal variation of temperature
const double DeltaTemp_Const[5] = { 0, 7.00, 11.00, 15.00, 14.50 };

// Seasonal variation of water vapor pressure
const double DeltaWaterVapPressure_Const[5] = { 0, 8.85, 7.24, 5.36, 3.39 };

// Seasonal variation of temperature lapse rate
const double DeltaTempLapsRt_Const[5] = { 0, 0.25e-3, 0.32e-3, 0.81e-3, 0.62e-3 };

// Seasonal variation of water vapor pressure lapse rate
const double DeltaWaterVapLaps_Const[5] = { 0, 0.33, 0.46, 0.74, 0.30 };


// Constant used in RTCA98 model
const double k1_Const = 77.604;
const double k2_Const = 382000;
const double Rd_Const = 287.054;
const double Gm_Const = 9.784;
const double g_Const = 9.80665;
const double RtcaMap1Pt001_Const = 1.001;
const double RtcaMapPt002001_Const = 0.002001;
const double TropoVerError_Const = 0.12;


// Constants used in both RTCA98 and Hopfield model
const double RSpecificAir_Const = 287.054;
const double Gravity_Const = 9.80665; // m/s2
const double Tropo9Pt784Const = 9.784;
const double RadianToDegree_Const = 180.0 / Pi_Const;
const double DegreetoRadian_Const = Pi_Const / 180.0;

// Temperature lapse rate for hopfield model
// and Sea level pressure(N/m^2), temperature(K) and density(kg/m^3)
const double aNot_Const = -6.5e-3;
const double aTwo_Const = 3e-3;
const double SeaLvlTemp_Const = 288.16;
const double SeaLvlPre_Const = 1.01325e5;
const double SeaLvlDensity_Const = 1.225;

// Geopotential height constants for standard atmosphere
const double hOne_Const = 1.1e+004;
const double hTwo_Const = 2.5e+004;
const double hThree_Const = 30000;

