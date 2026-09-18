// SIDD_1R.h: interface for the SIDD_1R class.
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

#ifndef SIDD_1R_H_
#define SIDD_1R_H_

#include "SIDD_Base.h"
#include "stat_1R.h"

using namespace std;

typedef	list<stat_1R> LstState;  

class SIDD_1R : public SIDD_Base  
{
protected:
	LstState Lst_OBE[MaxInitialWindowSize+1][3]; // lists of open base energy
	bool flag_minE1; // flagging if new min_E found in one run
	double min_RE[3]; // store minimum run energy (a + NI)
    long StatesInOneRun;

public:
	SIDD_1R();
	virtual ~SIDD_1R();
	void gen_OpenBaseEnergy();
	bool Search_Low1RE();
	
	void Update_Low1RE();
	void Show_Low1RE();
	double runs[3];
    double Prob[3];
    double exp_one[3];
    double exp_two[3][3];
    double exp_three[3][3][3];
    double exp_four[3][3][3][3];

	double exp_type_B;

    // V: This Prob_compete parameter was added for the competition branch
    // t=0 for melting, t=1 for Z-DNA, and t=2 for cruciforms
    double Prob_compete[3][3]; // V: This should be the probability of intersection?
    //  When [0][1], means you have runs of t1=0 and t2=1, melting and Z-DNA
    //  When [0][2], means you have t1=0 and t2=2, melting and cruciforms
    //  When [1][2], means you have t1=1 and t2=2, Z-DNA and cruciforms.
    //  Cases when [t][t] means having multiple runs of the same transition

    // V: Also these arrays were added for the competition branch
    // t=0 for melting, t=1 for Z-DNA, and t=2 for cruciforms
    double Prob_conditional[3][3]; // V: This should be the conditional probability
    double Prob_conditional_not[3][3]; // V: This should be the conditional probability of having [A] but not [B]


};

#endif // !defined(AFX_SIDD_1R_H__2688B961_9B45_44D5_8F6F_F45636E04460__INCLUDED_)
