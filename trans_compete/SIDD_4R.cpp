// SIDD_4R.cpp: implementation of the SIDD_4R class.
// This class is designed for four run derived from SIDD_Base
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

#ifndef SIDD_4R_CPP
#define SIDD_4R_CPP

#include "SIDD_4R.h"

// V: These new libraries were added
#include <vector>
#include <algorithm>
#include <cmath>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SIDD_4R::SIDD_4R():SIDD_3R()
{
	
	flag_minE4 = false;
	StatesInFourRuns = 0;
	sum_4RB = sum_4RG = 0.0;
}

SIDD_4R::~SIDD_4R()
{	
}

// search all lowest states for four runs
bool SIDD_4R::Search_Low4RE()
{
	double count = 0.0;
    	int par[3];
	for (int t1=0; t1<=Ns; t1++) {
        for (int t2=0; t2<=Ns; t2++) {
            for (int t3=0; t3<=Ns; t3++) {
                for (int t4=0; t4<=Ns; t4++) {
                    exp_four[t1][t2][t3][t4] = 0; 
                }
            }
        }
    }
	LstState::iterator i;
	LstState::iterator k;
	LstState::iterator j;
	LstState::iterator l;
	flag_minE4 = false;
	sum_4RB = sum_4RG = 0.0;
	int window_size = 200;
    double Epsilon = 0.0;
	for (int t1=0; t1<=2; t1++) {
		for (int t2=0; t2<=2; t2++) {
			for (int t3=0; t3<=2; t3++) {
				for (int t4=0; t4<=1; t4++) {					
					for(int n1 = 1; n1 <= window_size/4; n1++){
						//start here
							if(Lst_OBE[n1][t1].begin()->get_energy() + min_RE[t2] + min_RE[t3]+ min_RE[t4] + minGres >= max_E)
								continue;
						for(int n2 = n1; n2 <= (window_size - n1)/3 ; n2++){
                            if (n1==n2 && t1==1 && t2==0) 
                                continue;
								if(Lst_OBE[n1][t1].begin()->get_energy() + Lst_OBE[n2][t2].begin()->get_energy()+ min_RE[t3] + min_RE[t4] + minGres >= max_E)	
									continue;
							for(int n3 = n2; n3 <= (window_size - n1 - n2)/2; n3++){
								if (n1==n3 && t1==1 && t3==0) continue;
								if (n2==n3 && t2==1 && t3==0) continue;
									if(Lst_OBE[n1][t1].begin()->get_energy() + Lst_OBE[n2][t2].begin()->get_energy() + Lst_OBE[n3][t3].begin()->get_energy() + min_RE[t4] + minGres >= max_E)	
										continue;
								for(int n4 = n3; n4 <= (window_size - n1 - n2-n3); n4++){
									if (n1==n4 && t1==1 && t4==0) continue;
									if (n2==n4 && t2==1 && t4==0) continue;
									if (n3==n4 && t3==1 && t4==0) continue;
                                    for (int p=0; p<=2; p++)
                                        par[p]=n1*delta_fnc(t1,p)+n2*delta_fnc(t2,p)+n3*delta_fnc(t3,p)+n4*delta_fnc(t4,p);
                                    double Gr=calc_Gres(par[0],par[1],par[2],t1*delta_fnc(t1,1)+t2*delta_fnc(t2,1)+t3*delta_fnc(t3,1)+t4*delta_fnc(t4,1));
									if(Lst_OBE[n1][t1].begin()->get_energy() +Lst_OBE[n2][t2].begin()->get_energy() + Lst_OBE[n3][t3].begin()->get_energy() + Gr>= max_E)
										continue;
									for(i = Lst_OBE[n1][t1].begin(); i != Lst_OBE[n1][t1].end(); ++i){
										int p1 = i->get_pos1();
										double e1 = i->get_energy();
                                        if (e1 >= 10000)
                                            break;
										j = Lst_OBE[n2][t2].begin();
										if(n1 == n2 && t1==t2) ++j; // no repeat of the same group
										if(e1 + j->get_energy() + Lst_OBE[n3][t3].begin()->get_energy() + Lst_OBE[n4][t4].begin()->get_energy()+ Gr>= max_E)
											break;
										for(; j != Lst_OBE[n2][t3].end(); ++j){
											int p2 = j->get_pos1();
											double e2 = j->get_energy();
                                            if (e2 >= 10000) break;
											k = Lst_OBE[n3][t3].begin();
											if(n2 == n3 && t2==t3) ++k;
											if(e1 + e2 + k->get_energy() + Lst_OBE[n4][t4].begin()->get_energy()+Gr >= max_E)
												break;
											for(; k != Lst_OBE[n3][t3].end(); ++k){
												int p3 = k->get_pos1();
												double e3 = k->get_energy();
                                                if (e3 >= 10000) break;
												l = Lst_OBE[n4][t4].begin();
												if(n3 == n4 && t3==t4) ++l;
												if(e1 + e2 + e3 + l->get_energy() +Gr  >= max_E)
													break;
												for(; l != Lst_OBE[n4][t4].end(); ++l){
													int p4 = l->get_pos1();
													double e4 = l->get_energy();
                                                    if (e4 >= 10000)    break;
													double e = e1 + e2 + e3 + e4 +Gr;
													if(e >= max_E) break;
													if(!overlap(p1, p2, n1, n2) && !overlap(p1, p3, n1, n3) && !overlap(p1, p4, n1, n4) && !overlap(p2, p3, n2, n3) && !overlap(p2, p4, n2, n4) && !overlap(p3, p4, n3, n4)){
														if(min_E - e >=Epsilon){
															min_E = e;
															max_E = e + theta;
														}
														count++;
														
														if(write_profile){
															//double exponent = exp(-e/RT);
                                                            double exponent= exp(-e/RT + scaling_factor);  // V: scaling factor added

                                                            exp_four[t1][t2][t3][t4]+=exponent;
                                                            for (int t=0; t<=Ns; t++) 
                                                                runs[t] +=exponent*(delta_fnc(t,t1)+delta_fnc(t,t2)+delta_fnc(t,t3)+delta_fnc(t,t4));
															update_promatrix(p1, n1, t1, e, exponent);
															update_promatrix(p2, n2, t2, e, exponent);
															update_promatrix(p3, n3, t3, e, exponent);
															update_promatrix(p4, n4, t4, e, exponent);
															sum_4RB += exponent;
															sum_4RG = sum_4RG + e*exponent;

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
            }
        }
	}
	
	StatesInFourRuns = (long)count;
	ZsumB += sum_4RB;
	ZsumG += sum_4RG;
    totalFreq = StatesInOneRun + StatesInTwoRuns + StatesInThreeRuns+StatesInFourRuns;
    if (results) {
        cout << "Number of four-run states = " << count << endl;
        cout << "Total number of states = " << totalFreq << endl;

        // cout << "Total Partition Function without exp= " << ZsumB << endl;
        // V: The scaling factor needs to be removed from the partition function
        cout << "Total Partition Function = " << exp(-scaling_factor)*ZsumB << endl;
        cout << "Scaling factor " << scaling_factor << endl;

        cout << "Number of M_runs = " << runs[0]/ZsumB << endl;
        cout<< "Number of Z_runs = " << runs[1]/ZsumB << endl;
        cout<< "Number of C_runs = " << runs[2]/ZsumB << endl;
        double runs= (sum_1RB+2*sum_2RB+3*sum_3RB+4*sum_4RB)/ZsumB;
        cout << "Average number of runs = " << runs << endl;
    }

    // V: Here, they calculate the probabilities of having a particular transition.
    for (int t=0; t<=Ns; t++) {
        for (int t1=0; t1<=Ns; t1++) {
            for (int t2=0; t2<=Ns; t2++) {
                for (int t3=0; t3<=Ns; t3++) {
                    for (int t4=0; t4<=Ns; t4++) {
                        if (t1==t || t2==t || t3==t || t4==t)
                            Prob[t]+=exp_four[t1][t2][t3][t4]; 
                    }
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
                        for (int t4=0; t4<=Ns; t4++) {

                            // Let's build vectors of type of transitions, then check if these transitions are included
                            std::vector<int> firstVec = { t1, t2, t3, t4};
                            std::vector<int> secondVec = { ti, tj};

                            // Sort first vector
                            std::sort(firstVec.begin(), firstVec.end());
                            // Sort second vector
                            std::sort(secondVec.begin(), secondVec.end());

                            // Check if  all elements of a second vector exists in first vector
                            bool my_condition = std::includes(firstVec.begin(), firstVec.end(), secondVec.begin(), secondVec.end());

                            if (my_condition) {
                                Prob_compete[ti][tj] += exp_four[t1][t2][t3][t4];
                            }
                        }
                    }
                }
            }
        }
    }

    // V: Now, let's calculate the conditionals
    for (int t1=0; t1<=Ns; t1++) {
        for (int t2=0; t2<=Ns; t2++) {
            // V: The diagonal for Prob_conditional_not should be 0, right?
            Prob_conditional[t1][t2] = Prob_compete[t1][t2]/Prob[t2];
            // Probability of t1 given that we don't have t2
            double aa = Prob[t1]/ZsumB;
            double bb = Prob_compete[t1][t2]/ZsumB;
            double cc = Prob[t2]/ZsumB;

            //if (cc >= 1.0) Prob_conditional_not[t1][t2] = std::nan("NaN");
            // If the denominator is 0, then let's make the conditional probability 0.
            if (cc >= 1.0) Prob_conditional_not[t1][t2] = 0.0;
            else Prob_conditional_not[t1][t2] = (aa-bb)/(1.0-cc);  // If it isn't 0, then calculate it.
        }
    }

    if (results) {
        cout << "Prob_M = " <<  Prob[0]/ZsumB  << endl;
        cout<< "Prob_Z = "  << Prob[1]/ZsumB   << endl;
        cout << "Prob_C = " << Prob[2]/ZsumB  << endl;
    }

	flag_minE4 = 0;
	return flag_minE4;
}


#endif
