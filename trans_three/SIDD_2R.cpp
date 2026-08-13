// SIDD_2R.cpp: implementation of the SIDD_2R class.
// This class is designed for two runs derived from SIDD_1R
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
	long count = 0;
	LstState::iterator i;
	LstState::iterator j;
	flag_minE2 = false; // default value	
	sum_2RB = sum_2RG = 0.0;

    for(int n1 = MinWindowSize; n1 <= MaxWindowSize/2; n1++){
		if(Lst_OBE[n1].size() < 1)
			continue;
		if(Lst_OBE[n1].begin()->get_energy() + min_RE + minGres >= max_E)  
			continue;
		for(int n2 = n1; n2 <= MaxWindowSize - n1; n2++){
			if(Lst_OBE[n2].size() < 1)
				continue;
			if(Lst_OBE[n1].begin()->get_energy() + Lst_OBE[n2].begin()->get_energy() + Gres[n1+n2][1] >= max_E)
				continue;

            // V: I added this line to calculate the residual superhelical density
            //double a_residual = calc_alpha_res(n1+n2);

			for(i = Lst_OBE[n1].begin(); i != Lst_OBE[n1].end(); ++i){
				int p1 = i->get_pos1();	
				double e1 = i->get_energy();

                if (e1 >= 10000) break;				
				j = Lst_OBE[n2].begin();
				if(n1 == n2) ++j; // no repeat of the same group
				if(e1 + j->get_energy() + Gres[n1+n2][1]>= max_E) break;
				for(; j != Lst_OBE[n2].end(); ++j){							
					int p2 = j->get_pos1();
					double e2 = j->get_energy();
                    if (e2 >= 10000) break;
					double e = e1 + e2 + Gres[n1+n2][1]; //total energy for two runs
					if(e >= max_E) break;
					if(!overlap(p1, p2, n1, n2)){						
						if(e < min_E){
							min_E = e;
							max_E = min_E + theta;
							flag_minE2 = true;
						}
						if(e < max_E){														
                            if(write_profile){
                                double exponent= exp(-e/RT + scaling_factor);  // V: I added this
                                // double exponent= exp(-e/RT);
							 	update_promatrix(p1, n1, e, exponent);
								update_promatrix(p2, n2, e, exponent);
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
	StatesInTwoRuns = count;
	ZsumB += sum_2RB;
	ZsumG +=sum_2RG;
	
	if (results)
        cout << "Number of two-run states = " << count << endl;
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
