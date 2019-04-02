#pragma once
#include "c_main.h"

class CIrnss : public CMain
{
public:
	CIrnss();
	~CIrnss();
	void IrnssSat(struct TimeVar m_TVar, struct UserMot m_UserMotion, struct SatData *mp_IrnssSatData, double *mp_IrnssTransTime);

	struct UserMot m_IrnssUserMotion;
	struct TimeVar m_IrnssTVar;
	SatData *mp_IrnssSat;

	struct Sat
	{
		double  x_unit,
				y_unit,
				z_unit;
	};
	Sat *mp_Irnss; // Allocating memory space to pointer

private:
	int  IrnssYumaRead();
	void IrnssPosiVel(int sat_id);
	void IrnssRangeAndRate(int sat_id);
	void IrnssAziEle(int sat_id);
	int m_NoOfIrnssSat;
	static bool m_flag_irnss;

	struct IrnssYumaData
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
		af0s,
		af1s,
		week,
		semi_maj,
		m_ini;
	};

	IrnssYumaData *mp_IrnssYuma; 

	
	double  m_TimeDiff,
			m_TInt,
			*mp_TransitTime,
			m_SigTransmitTime;
	
};

