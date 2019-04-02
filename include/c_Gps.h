#pragma once
#include "c_main.h"

class CGps : public CMain
{
public:
	CGps();
	~CGps();
	void GpsSat(struct TimeVar m_TVar, struct UserMot m_UserMotion, struct SatData *mp_GpsSatData, double *mp_TransTime);

	struct UserMot m_GpsUserMotion;
	struct TimeVar m_GpsTVar;
	SatData *mp_GpsSat;

	struct Sat
	{
		double  x_unit,
				y_unit,
				z_unit;
	};
	Sat *mp_Gps; 
	
private:
	int  YumaRead();
	void PosiVel(int iter);
	void RangeAndRate(int iter);
	void AziEle(int iter);
	static bool m_flag_gps;
	int m_NoOfSat;

	struct GpsYumaData
	{
		double  sat_id,
		health,
		eccent,
		time,
		a_o_i,
		rt_r_asc,
		sqrt_semi_maj,
		r_asc,
		arg_of_pge,
		m_ini,
		af0s,
		af1s,
		week,
		semi_maj;
	};
	GpsYumaData *mp_Yuma;

	double  m_TimeDiff,
			m_TInt,
			*mp_TransitTime,
			m_SigTransmitTime;


};