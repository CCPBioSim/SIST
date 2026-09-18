// stat_1R.cpp: implementation of the stat_1R class.
// This class is for storing position of base pair and it's energy
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

#ifndef STAT_1R_CPP
#define STAT_1R_CPP

#include "stat_1R.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

stat_1R::stat_1R(int s, double e, double RT) //constructor
{
	start_pos1 = s;  //starting position
	energy = e;
}

stat_1R::~stat_1R()  //destructor ???
{

}



bool stat_1R::operator<(const stat_1R& x1)const
{

	return (energy < x1.get_energy());

}

bool stat_1R::operator==(const stat_1R& x1)const
{

	return (energy == x1.get_energy());

}

bool stat_1R::operator>(const stat_1R& x1)const
{

	return (energy > x1.get_energy());

}

#endif
