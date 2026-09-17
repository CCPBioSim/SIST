// stat_1R.h: interface for the stat_1R class.
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

#include <math.h>
#ifndef STAT_1R_H_
#define STAT_1R_H_

class stat_1R  
{
protected:
	int start_pos1;
	double energy;
public:
	stat_1R(int start_position, double energy, double rt);
	virtual ~stat_1R();
	int get_pos1() const{return start_pos1;}
	bool operator <(const stat_1R&) const;
	bool operator ==(const stat_1R&)const;
	bool operator >(const stat_1R&) const;
	double get_energy()const {return energy;}
};

#endif // !defined(AFX_STAT_1R_H__ACD4635F_79AA_463A_A66C_6B1CA4476EE0__INCLUDED_)
