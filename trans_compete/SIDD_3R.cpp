// SIDD_3R.cpp: implementation of the SIDD_3R class.
// This class is designed for three run derived from SIDD_Base
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

#ifndef SIDD_3R_CPP
#define SIDD_3R_CPP

#include "SIDD_3R.h"

#include <vector>
#include <algorithm>


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SIDD_3R::SIDD_3R():SIDD_2R()
{
	
	flag_minE3 = false;
	StatesInThreeRuns = 0;
	sum_3RB = sum_3RG = 0.0;
}

SIDD_3R::~SIDD_3R()
{	
}

// search all lowest states for three runs
bool SIDD_3R::Search_Low3RE()
{
	
	double count = 0.0;
    	int par[3];
	for (int t1=0; t1<=Ns; t1++) {
        for (int t2=0; t2<=Ns; t2++) {
            for (int t3=0; t3<=Ns; t3++) {
            exp_three[t1][t2][t3] = 0; 
            }
        }
    }
	LstState::iterator i;
	LstState::iterator k;
	LstState::iterator j;
	flag_minE3 = false;
	sum_3RB = sum_3RG = 0.0;
	int window_size = 200;
    double Epsilon = 0.0;
	for (int t1=0; t1<=2; t1++) {
		for (int t2=0; t2<=2; t2++) {
			for (int t3=0; t3<=2; t3++) {				
				for(int n1 = 1; n1 <= window_size/3; n1++){
					if(Lst_OBE[n1][t1].begin()->get_energy() + min_RE[t2]+min_RE[t3] + minGres >= max_E)
						continue;
					for(int n2 = n1; n2 <= (window_size - n1)/2 ; n2++){
						if (n1==n2 && t1==1 && t2==0) continue;
						if(Lst_OBE[n1][t1].begin()->get_energy() + Lst_OBE[n2][t2].begin()->get_energy()+ min_RE[t3] + minGres >= max_E)
							continue;
						for(int n3 = n2; n3 <= window_size - n1 - n2; n3++){
							if (n1==n3 && t1==1 && t3==0) continue;
							if (n2==n3 && t2==1 && t3==0) continue;
                            for (int p=0; p<=2; p++)
                                par[p]=n1*delta_fnc(t1,p)+n2*delta_fnc(t2,p)+n3*delta_fnc(t3,p);
                            double Gr=calc_Gres(par[0],par[1],par[2],t1*delta_fnc(t1,1)+t2*delta_fnc(t2,1)+t3*delta_fnc(t3,1));
							if(Lst_OBE[n1][t1].begin()->get_energy() + Lst_OBE[n2][t2].begin()->get_energy() + Lst_OBE[n3][t3].begin()->get_energy() + Gr>= max_E)
								continue;
							for(i = Lst_OBE[n1][t1].begin(); i != Lst_OBE[n1][t1].end(); ++i){
								int p1 = i->get_pos1();
								double e1 = i->get_energy();
                                if (e1 >= 10000) break;								
								j = Lst_OBE[n2][t2].begin();
								if(n1 == n2 && t1==t2) ++j; // no repeat of the same group
								if(e1 + j->get_energy() + Lst_OBE[n3][t3].begin()->get_energy()+Gr >= max_E)
									break;
								for(; j != Lst_OBE[n2][t2].end(); ++j){
									int p2 = j->get_pos1();
									double e2 = j->get_energy();
                                    if (e2 >= 10000) break;
									k = Lst_OBE[n3][t3].begin();
									if(n2 == n3 && t2==t3) ++k;
									if(e1 + e2 + k->get_energy() + Gr >= max_E)
										break;
									for(; k != Lst_OBE[n3][t3].end(); ++k){
										int p3 = k->get_pos1();
										double e3 = k->get_energy();
                                        if (e3 >= 10000)    break;
										double e = e1 + e2 + e3 + Gr;
										if(e >= max_E) break;
										if(!overlap(p1, p2, n1, n2) && !overlap(p1, p3, n1, n3) && !overlap(p2, p3, n2, n3)){
											if(min_E - e >=Epsilon){
												min_E = e;
												max_E = e + theta;
											}
											count++;
											if(write_profile){
												// double exponent = exp(-e/RT);
                                                double exponent= exp(-e/RT + scaling_factor);  // V: scaling factor added
												exp_three[t1][t2][t3] +=exponent;
                                                for (int t=0; t<=Ns; t++) 
                                                    runs[t] +=exponent*(delta_fnc(t,t1)+delta_fnc(t,t2)+delta_fnc(t,t3));
												update_promatrix(p1, n1, t1, e, exponent);
												update_promatrix(p2, n2, t2, e, exponent);
												update_promatrix(p3, n3, t3, e, exponent);
												sum_3RB += exponent;
												sum_3RG = sum_3RG + e*exponent;
											}
										}
									}
								}
							}
						}
					}
				}
            }
        }
	}
	
	if (results)
        cout << "Number of three-run states = " << count << endl;
	StatesInThreeRuns = (long)count;
	ZsumB += sum_3RB;
	ZsumG += sum_3RG;
    for (int t=0; t<=Ns; t++) {
        for (int t1=0; t1<=Ns; t1++) {
            for (int t2=0; t2<=Ns; t2++) {
                for (int t3=0; t3<=Ns; t3++) {
                    if (t1==t || t2==t || t3==t)
                        Prob[t]+=exp_three[t1][t2][t3]; 
                }
            }
        }
    }
    // V: This was added for the competition of transitions:
    for (int ti=0; ti<=Ns; ti++) {
        for (int tj=0; tj<=Ns; tj++) {
            for (int t1=0; t1<=Ns; t1++) {
                for (int t2=0; t2<=Ns; t2++) {
                    for (int t3=0; t3<=Ns; t3++) {

                        // Let's build vectors of type of transitions, then check if these transitions are included
                        std::vector<int> firstVec = { t1, t2, t3};
                        std::vector<int> secondVec = { ti, tj};

                        // Sort first vector
                        std::sort(firstVec.begin(), firstVec.end());
                        // Sort second vector
                        std::sort(secondVec.begin(), secondVec.end());

                        // Check if  all elements of a second vector exists in first vector
                        bool my_condition = std::includes(firstVec.begin(), firstVec.end(), secondVec.begin(), secondVec.end());

                        if (my_condition) {
                            Prob_compete[ti][tj] += exp_three[t1][t2][t3];
                            //cout << " ti "<<ti<< " tj "<<tj<< " t1 "<<t1<< " t2 "<<t2<< " t3 "<<t3<<endl;
                        }
                    }
                }
            }
        }
    }
    flag_minE3 = 0;
return flag_minE3;
}

#endif
