// SIDD_3R.cpp: implementation of the SIDD_3R class.
// This class is designed for three runs derived from SIDD_2R
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
	LstState::iterator i;
	LstState::iterator k;
	LstState::iterator j;
	flag_minE3 = false;
	sum_3RB = sum_3RG = 0.0;

    int window_size = 200;
    double Epsilon = 0.0;
	for(int n1 = 1; n1 <= window_size/3; n1++){
		if(Lst_OBE[n1].begin()->get_energy() + 2*min_RE + minGres >= max_E) continue;
		for(int n2 = n1; n2 <= (window_size - n1)/2 ; n2++){
			if(Lst_OBE[n1].begin()->get_energy() + Lst_OBE[n2].begin()->get_energy()+ min_RE + minGres >= max_E) continue;
			for(int n3 = n2; n3 <= window_size - n1 - n2; n3++){
				if(Lst_OBE[n1].begin()->get_energy() + Lst_OBE[n2].begin()->get_energy() + Lst_OBE[n3].begin()->get_energy() + Gres[n1+n2+n3][2]>= max_E) continue;

                // V: I added this line to calculate the residual superhelical density
                //double a_residual = calc_alpha_res(n1+n2+n3);

                for(i = Lst_OBE[n1].begin(); i != Lst_OBE[n1].end(); ++i){
					int p1 = i->get_pos1();
					double e1 = i->get_energy();
                    if (e1 >= 10000) break;
					j = Lst_OBE[n2].begin();
					if(n1 == n2) ++j; // no repeat of the same group
					if(e1 + j->get_energy() + Lst_OBE[n3].begin()->get_energy() + Gres[n1+n2+n3][2] >= max_E) break;
					for(; j != Lst_OBE[n2].end(); ++j){
						int p2 = j->get_pos1();
						double e2 = j->get_energy();
                        if (e2 >= 10000) break;
						k = Lst_OBE[n3].begin();
						if(n2 == n3) ++k;
						if(e1 + e2 + k->get_energy() + Gres[n1+n2+n3][2] >= max_E) break;
						for(; k != Lst_OBE[n3].end(); ++k){
							int p3 = k->get_pos1();
							double e3 = k->get_energy();
                            if (e3 >= 10000) break;
							double e = e1 + e2 + e3 +Gres[n1+n2+n3][2];
							if(e >= max_E) break;
                            //checking for overlaps
							if(!overlap(p1, p2, n1, n2) && !overlap(p1, p3, n1, n3) && !overlap(p2, p3, n2, n3)){
								if(min_E - e >=Epsilon){
									min_E = e;
									max_E = e + theta;
								}
								count++;
								if(write_profile){
									// double exponent = exp(-e/RT);
                                    double exponent= exp(-e/RT + scaling_factor);  // V: I added this

                                    update_promatrix(p1, n1, e, exponent);
									update_promatrix(p2, n2, e, exponent);
									update_promatrix(p3, n3, e, exponent);
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
	if (results)
        cout << "Number of three-run states = " << count << endl;
	StatesInThreeRuns = (long)count;
	ZsumB += sum_3RB;
	ZsumG += sum_3RG;
    return flag_minE3;
}

#endif
