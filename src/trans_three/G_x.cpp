// G_x.cpp: implementation of the G_x class.
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



#ifndef G_x_CPP
#define G_x_CPP

#include "G_x.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

G_x::G_x()
{
	sum_xG = 0.0;
	sum_xB = 0.0; 
	ave_Gx = 0.0;
	px = 0.0;
}

void G_x::reset()
{
	sum_xG = 0.0; 
	sum_xB = 0.0;
	ave_Gx = 0.0;
	px = 0.0;
}


void G_x::add(double gs, double exponent, double rt, double lastG, double lastB)
{
        if(gs == -1.0 && rt == -1.0){
                sum_xG = sum_xG + lastG;
                sum_xB = sum_xB + lastB;
        }
        else{
                sum_xG = sum_xG + gs*exponent + lastG;
                sum_xB = sum_xB + exponent + lastB;
        }
}


// HERE IS THE DIVISION THAT MAKES IT P=1.0, when ZsumB == 0.0. We need another way to figure out this.
void G_x::calc(double ave_Gs, double ZsumB) //input: average energy and partition function
{
	if (ZsumB == 0.0)
		px = 1.0;
	else
		px = sum_xB / ZsumB;  
	
	if(sum_xB == 0.0)
		ave_Gx = INFINITE_Gx;
	else
		ave_Gx = sum_xG / sum_xB - ave_Gs;
}

G_x::~G_x()
{

}

#endif
