// SIDD_4R.h: interface for the SIDD_4R class.
// This class is designed for four runs derived from SIDD_3R
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


#ifndef SIDD_4R_H_
#define SIDD_4R_H_

#include "SIDD_3R.h"

class SIDD_4R : public SIDD_3R  
{
private:
	bool flag_minE4;
	long StatesInFourRuns;

public:
	SIDD_4R();
	virtual ~SIDD_4R();
	bool Search_Low4RE();
};

#endif // !defined(AFX_SIDD_3R_H__B5CD0095_10DE_4C22_A817_256ABBAD9425__INCLUDED_)
