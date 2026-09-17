// SIDD_2R.h: interface for the SIDD_2R class.
//
// program to implement algorithm developed by Craig Benham
//
// author: Chengpeng Bi
// modifiers: Dina Zhabinskaya, Sally Madden, Ian Korf
// compiler: g++
//
// This is a test version. The author has no responsibility for any outcome
// that incurs when you make any trial run.
//
// UC Davis Genome Center
//////////////////////////////////////////////////////////////////////////////

#ifndef SIDD_2R_H_
#define SIDD_2R_H_

#include "SIDD_1R.h"


class SIDD_2R : public SIDD_1R  
{
protected:
	
	bool flag_minE2;
	long StatesInTwoRuns;

public:
	SIDD_2R();
	virtual ~SIDD_2R();
	
	// two-run states
	bool overlap(int p1, int p2, int r1, int r2);
	bool Search_Low2RE();

};

#endif // !defined(AFX_SIDD_2R_H__79785266_E18C_42AD_8C80_6BC6C04A20DD__INCLUDED_)
