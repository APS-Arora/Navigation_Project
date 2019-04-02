/********************************************************************************************************
* Description            : Main program to call all the other classes
* Specific library calls : None
* Classes                : 1. c_class
* Assumptions            : None
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
#include "c_Gps.h"
#include "c_Irnss.h"
#include "GlobalFunction.h"
#include "ConstVar.h"
#include "c_main.h"

using namespace std;
using namespace Eigen;



int main()
{

	// Create object for c_main class
	CMain CMain_Obj;

	// Function call from c_main class
	CMain_Obj.MainFunc();
	

	
	return 0;
}
