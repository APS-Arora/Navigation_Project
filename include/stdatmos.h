#pragma once
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
using namespace std;
//#include <string.h>

// P H Y S I C A L   C O N S T A N T S
const double FT2METERS = 0.3048;     // mult. ft. to get meters (exact)
const double KELVIN2RANKINE = 1.8;              // mult kelvins to get deg R
const double PSF2NSM = 47.880258;      // mult lb/sq.ft to get N/sq.m
const double SCF2KCM = 515.379;    // mult slugs/cu.ft to get kg/cu.m
const double TZERO = 288.15;      // sea level temperature, kelvins
const double PZERO = 101325.0;        // sea-level pressure, N/sq.m
const double RHOZERO = 1.225;           // sea level density, kg/cu.m
const double AZERO = 340.294;    // sea-level speed of sound, m/sec

// F U N C T I O N   P R O T O T Y P E S
void SimpleAtmosphere(const double, double&, double&, double&);
void Atmosphere(const double, double&, double&, double&);
double MetricViscosity(const double);