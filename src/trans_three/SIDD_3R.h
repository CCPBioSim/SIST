// SIDD_3R.h: interface for the SIDD_3R class.
// This class is designed for three runs derived from SIDD_2R
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

#ifndef SIDD_3R_H_
#define SIDD_3R_H_

#include "SIDD_2R.h"

class SIDD_3R : public SIDD_2R  
{
protected:
	bool flag_minE3;
	long StatesInThreeRuns;

public:
	SIDD_3R();
	virtual ~SIDD_3R();
	bool Search_Low3RE();
};

#endif // !defined(AFX_SIDD_3R_H__B5CD0095_10DE_4C22_A817_256ABBAD9425__INCLUDED_)
