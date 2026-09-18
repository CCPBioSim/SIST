// SIDD_2R.cpp: implementation of the SIDD_2R class.
// This class is designed for two run derived from SIDD_Base
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

#ifndef SIDD_2R_CPP
#define SIDD_2R_CPP

#include "SIDD_2R.h"

#include <vector>
#include <algorithm>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SIDD_2R::SIDD_2R():SIDD_1R()
{
	flag_minE2 = false;
	StatesInTwoRuns = 0;
	sum_2RB = sum_2RG = 0.0;
	
}

SIDD_2R::~SIDD_2R()
{	
}

// search all lowest states for two runs
bool SIDD_2R::Search_Low2RE()
{	
    int par[3];
	long count = 0;
	for (int t1=0; t1<=Ns; t1++) {
        for (int t2=0; t2<=Ns; t2++) {
            exp_two[t1][t2] = 0;

            // V: This was added for the probability of having at least two runs of two transitions
            Prob_compete[t1][t2]=0;
            // V: These will be the conditionals
            Prob_conditional[t1][t2]=0;
            Prob_conditional_not[t1][t2]=0;

        }
    }
	LstState::iterator i;
	LstState::iterator j;
	flag_minE2 = false; // default value

	sum_2RB = sum_2RG = 0.0;
	for (int t1=0; t1<=Ns; t1++) {
		for (int t2=0; t2<=Ns; t2++) {
			for(int n1 = MinWindowSize; n1 <= MaxWindowSize/2; n1++){
				if(Lst_OBE[n1][t1].size() < 1) 
					continue;
				if(Lst_OBE[n1][t1].begin()->get_energy() + min_RE[t2] + minGres >= max_E)  
					continue;
				for(int n2 = n1; n2 <= MaxWindowSize - n1; n2++){
					if (n1==n2 && t1==1 && t2==0) continue;
                    if (n1==n2 && t1==2 && t2==0) continue;
                    if (n1==n2 && t1==2 && t2==1) continue;
					if(Lst_OBE[n2][t2].size() < 1) continue;
                    for (int p=0; p<=2; p++)
                       par[p]=n1*delta_fnc(t1,p)+n2*delta_fnc(t2,p);
                    double Gr=calc_Gres(par[0],par[1],par[2],t1*delta_fnc(t1,1)+t2*delta_fnc(t2,1));
                    if(Lst_OBE[n1][t1].begin()->get_energy() + Lst_OBE[n2][t2].begin()->get_energy() + Gr >= max_E)
						continue;
					for(i = Lst_OBE[n1][t1].begin(); i != Lst_OBE[n1][t1].end(); ++i){		
						int p1 = i->get_pos1();	
						double e1 = i->get_energy();
                        if (e1 >= 10000) break;
						j = Lst_OBE[n2][t2].begin();
						if(n1 == n2 && t1==t2) ++j; // no repeat of the same group
						if (e1 + j->get_energy() + Gr >= max_E)  
							break;			
						for(; j != Lst_OBE[n2][t2].end(); ++j){		
							int p2 = j->get_pos1();
							double e2 = j->get_energy();
                            if (e2 >= 10000) break;
							double e = e1 + e2 + Gr; //total energy for two runs
							if (e >= max_E) break;
							if(!overlap(p1, p2, n1, n2)){						
								if(e < min_E){
									min_E = e;
									max_E = min_E + theta;
									flag_minE2 = true;
								}
								if(e < max_E){										
									if(write_profile){
                                        // double exponent= exp(-e/RT);
                                        // V: the scaling factor is added
                                        double exponent= exp(-e/RT + scaling_factor);

										exp_two[t1][t2] += exponent;
                                        for (int t=0; t<=Ns; t++) 
                                            runs[t] +=exponent*(delta_fnc(t,t1)+delta_fnc(t,t2));
										update_promatrix(p1, n1, t1, e, exponent);
										update_promatrix(p2, n2, t2, e, exponent);
										sum_2RB += exponent;
										sum_2RG = sum_2RG + e*exponent;
									}
									count++;

                                }
							}
						}
					}
				}
			}
		}
	}
	
	StatesInTwoRuns = count;
	ZsumB += sum_2RB;
	ZsumG +=sum_2RG;
    if (results)
        cout << "Number of two-run states = " << count << endl;
    for (int t=0; t<=Ns; t++) {
        for (int t1=0; t1<=Ns; t1++) {
            for (int t2=0; t2<=Ns; t2++) {
                if (t1==t || t2==t)
                    Prob[t] += exp_two[t1][t2];
            }
        }
    }
    // V: This counts for states where you have two different structural transitions, ti and tj.
    // V: Watch out for the diagonal, it doesn't tell you the prob of having one run of another of the same transition. It doesn't even give you p(M)
    for (int ti=0; ti<=Ns; ti++) {
        for (int tj=0; tj<=Ns; tj++) {
            for (int t1=0; t1<=Ns; t1++) {
                for (int t2=0; t2<=Ns; t2++) {

                    // Let's build vectors of type of transitions, then check if these transitions are included
                    std::vector<int> firstVec = { t1, t2};
                    std::vector<int> secondVec = { ti, tj};

                    // Sort first vector
                    std::sort(firstVec.begin(), firstVec.end());
                    // Sort second vector
                    std::sort(secondVec.begin(), secondVec.end());

                    // Check if  all elements of a second vector exists in first vector
                    bool my_condition = std::includes(firstVec.begin(), firstVec.end(), secondVec.begin(), secondVec.end());

                    if (my_condition) {
                        Prob_compete[ti][tj] += exp_two[t1][t2];
                    }
                }
            }
        }
    }

    return flag_minE2;
}




// checking if two regions are overlapping
bool SIDD_2R::overlap(int p1, int p2, int r1, int r2)
{
	
	if(r1 + r2 > MaxWindowSize){ 
		cerr << "out of window size.\n";
		return false;
	}
	
	if(p2 + r2 <= length_seq - 1 && p1 + r1 <= length_seq - 1){
		if(p1 + r1 < p2 - 1 || p2 + r2 < p1 - 1) // at least one gap
			return false; // no overlapping
		else
			return true;
	}
	else if(p2 + r2 > length_seq - 1 && p1 + r1 <= length_seq - 1){
		if(p1 + r1 < p2 - 1 && r2 - (length_seq - 1 - p2) < p1 - 1)
			return false; // no overlapping
		else
			return true;
	}
	else if(p1 + r1 > length_seq - 1 && p2 + r2 <= length_seq - 1){
		if(p2 + r2 < p1 - 1 && r1 - (length_seq - 1 - p1) < p2 - 1)
			return false; // no overlapping
		else
			return true;
	}
	else
		return true; // overlapping
}

#endif
