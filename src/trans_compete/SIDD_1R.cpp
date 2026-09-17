// SIDD_1R.cpp: implementation of the SIDD_1R class.
// This class is designed for one run derived from SIDD_Base
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

#ifndef SIDD_1R_CPP
#define SIDD_1R_CPP

#include "SIDD_1R.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SIDD_1R::SIDD_1R():SIDD_Base() {
	for (int t=0; t<=Ns; t++) 
		min_RE[t] = min_E;  //minimun energy from the continuous part 
        flag_minE1 = false;
        StatesInOneRun = 0;
}
SIDD_1R::~SIDD_1R()
{
}

void SIDD_1R::gen_OpenBaseEnergy()
{
    double e;
    //t=0 for melting, t=1 for Z-DNA, and t=2 for cruciforms
	for (int t=0; t<=Ns; t++) { 
		for(int i = MinWindowSize; i <= MaxWindowSize; i++){ // window size: 1 - MaxWindowSizez
			for(int j = 0; j < length_seq; j++){ // start position to length_seq - 1
                // V: a - is the nucleation energy to initiate run of transition t
                // V: calc_OPenBasesEnergy. The function is in SIDD_base. It calculates the energy associated
                //    with transition t given the run at position startp:j and with size n:i
                e = a[t] + calc_OPenBasesEnergy(j, i, t);
                // V: It creates object of type stat_1R. It seems it just stores the values of position j,
                //    energy e and RT
                stat_1R  s1r(j, e, RT);	 //input j=starting position, e=energy, RT=constant
                // V: It then appends this object the Lst_OBE array/list.
                //    This Lst_OBE seems to be a list made for each run/window of size i with transition t.
				Lst_OBE[i][t].push_back(s1r);
				if(e < min_RE[t]) min_RE[t] = e; 
			}
		}
	}
	
	// sorting each list of OBE
    // V: I see, for the list of windows/runs of size k of transition t=0-melting, t=1-Z, t=2 cruciform, it sorts
    //    from low to high. I think this is useful for filtering high energy states...
	for (int t=0; t<=Ns; t++) {
		for(int k = MinWindowSize; k <= MaxWindowSize; k++){	
			Lst_OBE[k][t].sort();
		}
	}
}


bool SIDD_1R::Search_Low1RE()
{
	long count = 0;
	flag_minE1 = false;
	sum_1RG = sum_1RB = 0.0;
	for (int t=0; t<=Ns; t++) {
		exp_one[t] = 0; 
		runs[t]=0; 
        Prob[t]=0;
    }

    LstState::iterator j;
	//collecting one run states withing the energy threshold
    // V: t is the index indicating transition, t=0 melting, t=1 Z-DNA and t=2 cruciform
	for (int t=0; t<=Ns; t++) {
		for(int i = MinWindowSize; i <= MaxWindowSize; i++){ // window size: 1 - MaxWindowSize
            // V: This iterates over the list which is composed of states with run of size i with transition t.
			for(j = Lst_OBE[i][t].begin(); j != Lst_OBE[i][t].end(); ++j){
                // V: This line retrieves the energy from the iterator j.
				double e = j->get_energy(); 
               if (e >= 10000)
                   break;
				int pos =j->get_pos1();
                // V: Calculates the superhelical energy? And adds it. I think it is the residual energy.
                //    I think it makes sense because each transition absorbs certain amount of superhelicity /
                //    free energy
				e += calc_Gres(i*delta_fnc(t,0),i*delta_fnc(t,1),i*delta_fnc(t,2),t*delta_fnc(t,1));
                //  V: This pretty much does what it says. it doesn't count high energy states
				if(e >= max_E) break; //outside of threshold
               				count++; //counting the number of one run states
                // V: If minimum energy found, register it.
				if(e < min_E){ 
					min_E = e; //update min_E
                    // V: theta = threshold. You always look for states that are within the minimum energy state
                    //    plus the threshold.
					max_E = e + theta;
					min_WS = i; // the window size corresponding to the minimum energy
					flag_minE1 = true;
				}
                // V: This is where it stores some values
                // V: So, basically what this does is that the exponent of the energy e retrieved from
                //    Lst_OBE (see above) is calculated.
                //    Then, these exponents are added to the quantitiy exp_one[t] according transition t.
                //    The number of runs[t] per transition t starts collecting this exponents as well.
                //    Then the matrix update_promatrix is updated - I still don't know what promatrix  is, but I think
                //    it is a vector list with coordinates promatrix[window_size][structural][position] and contains
                //    startp:position, n:window_size, struc:structural_transition, x:free_energy,
                //    bzfactor:Boltzman factor
                //    I believe each entry of promatrix has the form of G_x, but I'm not that sure
                //    sum_1RG sums all the states within one run.
                //    sum_1RG sums all Boltzman factors within one run
				if(write_profile){
                    // double exponent= exp(-e/RT);
                    double exponent= exp(-e/RT + scaling_factor);  // V: I added this scaling factor
                    // double exponent_test = exp(-e/RT + 200.0);
                    // V: exp_one[] is a vector of dimension 3, that sums all the exponents of runs 1.
                    //    it is used for calculating the probability of transition (see below).
                    //    that probability sums all these calculated exponents.
                    exp_one[t] += exponent;
                    // V: It doesn't make a difference if runs is inside the loop?
                    //    But I guess it counts the number of runs.
					runs[t]=exp_one[t];
					update_promatrix(pos, i, t, e, exponent);
                    sum_1RB += exponent;
					// sum_1RB += exponent_test;
					sum_1RG = sum_1RG + e*exponent; 
				}                  
			}
		}
	}
	
	StatesInOneRun = count;
	// adding the zero open bases term to the partition function
	double closed = exp(-alpha*alpha*K/2/RT + scaling_factor);
    // V: Note that the zero runs state is added to the partition function, but not to the probabilities.
	ZsumB = sum_1RB + closed;
    // V: Note that ZsumG has the form like sum_1RG (up there). Which is sum_1RG + e*exponent, where e = alpha*alpha*K/2
	ZsumG = sum_1RG+alpha*alpha*K/2*closed;

    if(write_profile && results) {
		cout << "Number of one-run states = " << count << endl;        
        for (int t=0; t<=Ns; t++) {
            Prob[t] += exp_one[t];
        }
    }
    return flag_minE1;
}

#endif

