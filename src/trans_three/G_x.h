// G_x.h: interface for the G_x class.
// This class is designed for storing profiles e.g. p(x), G(x), A(x)
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

#ifndef G_X_H_
#define G_X_H_

#include <math.h>

const double INFINITE_Gx = -10000.0;
class G_x
{
private:
	double sum_xG; // total energy with states of x opening
	double sum_xB; // sum of Boltzman's factors with states of x open
    double ave_Gx; // average free energy of opening base x
	double px;	   // probability of opening base x

public:
	G_x();
	void reset();
    void add(double gs, double exponent, double rt, double lastG = 0.0, double lastB = 0.0); // add one state
	void calc(double ave_Gs, double zsumb); // computing p(x) and G(x)
	double get_sum_xG(){return sum_xG;}
	double get_ave_Gx(){return ave_Gx;}
	void set_ave_Gx(double x){ave_Gx = x;}
	double get_sum_xB(){return sum_xB;}
	double get_px(){return px;}
	virtual ~G_x();

};

#endif 
