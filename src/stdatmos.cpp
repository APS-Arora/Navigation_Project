// PROGRAM stdatmos modified version of Tables whose description is below

// PROGRAM Tables;                                    { \atmtable\Tables.cpp
// --------------------------------------------------------------------------
// PURPOSE - Make tables of atmospheric properties
// AUTHOR  - Ralph L. Carmichael, Public Domain Aeronautical Software
// REFERENCE - U.S. Standard Atmosphere, 1976. U.S. Govt Printing Office
// REVISION HISTORY
//   DATE  VERS PERSON  STATEMENT OF CHANGES
// 26Feb95  1.0   RLC   Assembled several old codes
//  7Jul95  1.1   RLC   Added viscosity calculations
//  6Aug95  1.2   RLC   Replaced 1962 tables with 1976
// 25Aug95  1.3   RLC   Added MODIFIER
// 19Sep95  1.4   RLC   Added a little precision to pressure table
// 31Oct01  1.5   RLC   Renamed the output files
// 31Jul08  1.6   RLC   Adapted to standard library namespace
// --------------------------------------------------------------------------
#include "stdatmos.h"

// ==========================================================================
void Atmosphere(const double  alt,                  // geometric altitude, km.
		double& sigma,           // density/sea-level standard density
		double& delta,         // pressure/sea-level standard pressure
		double& theta)   // temperature/sea-level standard temperature
// Compute the temperature,density, and pressure in the standard atmosphere
// Correct to 86 km.  Only approximate thereafter.
{
  const double REARTH=6369.0;    // radius of the Earth (km)
  const double GMR = 34.163195;
  const int NTAB = 8;
  int i,j,k;

  static double htab[NTAB] = {0.0,  11.0, 20.0, 32.0, 47.0,
			     51.0, 71.0, 84.852 };
  static double ttab[NTAB] = { 288.15, 216.65, 216.65, 228.65, 270.65,
			      270.65, 214.65, 186.946 };
  static double ptab[NTAB] = { 1.0, 2.2336110E-1, 5.4032950E-2, 8.5666784E-3,
     1.0945601E-3, 6.6063531E-4, 3.9046834E-5, 3.68501E-6 };
  static double gtab[NTAB] = { -6.5, 0, 1.0, 2.8, 0, -2.8, -2.0, 0 };

  double h=alt*REARTH/(alt+REARTH);     //  geometric to geopotential altitude

  i=0; j=NTAB-1;  // starting values for binary search
  do
    {
      k=(i+j)/2;
      if (h < htab[k]) j=k; else i=k;
    }  while (j > i+1);

  double tgrad=gtab[i];                      // temp. gradient of local layer
  double tbase=ttab[i];                      // base temp. of local layer
  double deltah=h-htab[i];                   // height above local base
  double tlocal=tbase+tgrad*deltah;          // local temperature
  theta=tlocal/ttab[0];                                  // temperature ratio

  if (0.0 == tgrad)                                         // pressure ratio
    delta=ptab[i]*exp(-GMR*deltah/tbase);
  else
    delta=ptab[i]*pow(tbase/tlocal, GMR/tgrad);

  sigma=delta/theta;                                        //  density ratio
}   // ------------------------------------------- End of function Atmosphere

// ==========================================================================
void SimpleAtmosphere(
	const double alt,                           // geometric altitude, km.
	double& sigma,                   // density/sea-level standard density
	double& delta,                 // pressure/sea-level standard pressure
	double& theta)           // temperature/sea-level standard temperature

// Compute the temperature,density, and pressure in the standard atmosphere
// Correct to 20 km.  Only approximate thereafter.
{
  const double REARTH = 6369.0;    // radius of the Earth (km)
  const double GMR    = 34.163195;   // gas constant


  double h=alt*REARTH/(alt+REARTH);     //  geometric to geopotential altitude

  if (h<11.0)
    {                                                          // Troposphere
      theta=(288.15-6.5*h)/288.15;
      delta=pow(theta, GMR/6.5);
    }
  else
    {                                                         // Stratosphere
      theta=216.65/288.15;
      delta=0.2233611*exp(-GMR*(h-11.0)/216.65);
    }

  sigma=delta/theta;
}   // ------------------------------------- End of function SimpleAtmosphere

// ==========================================================================
double MetricViscosity(const double theta)
{
  const double TZERO = 288.15;               // sea level temperature, kelvins
  const double BETAVISC = 1.458E-6;     // constant, N-sec/(sq.m-sqrt(kelvins)
  const double SUTH = 110.4;                 // Sutherland's constant, kelvins

  double t=theta*TZERO;                              // temperature in kelvins
  return BETAVISC*sqrt(t*t*t)/(t+SUTH);            // viscosity in kg/(m-sec)
}   // -------------------------------------- End of function MetricViscosity
