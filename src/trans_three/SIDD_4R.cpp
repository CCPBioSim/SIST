// SIDD_4R.cpp: implementation of the SIDD_4R class.
// This class is designed for four runs derived from SIDD_3R
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

// search all lowest states for three runs
bool SIDD_4R::Search_Low4RE()
{
	double count = 0.0;
	LstState::iterator i;
	LstState::iterator k;
	LstState::iterator j;
	LstState::iterator l;
	flag_minE4 = false;
	sum_4RB = sum_4RG = 0.0;
	int window_size = 200;
    double Epsilon = 0.0;
	for(int n1 = 1; n1 <= window_size/4; n1++){
		if(Lst_OBE[n1].begin()->get_energy() + 3*min_RE + minGres >= max_E)
			continue;
		for(int n2 = n1; n2 <= (window_size - n1)/3 ; n2++){
			if(Lst_OBE[n1].begin()->get_energy() + Lst_OBE[n2].begin()->get_energy()+ 2*min_RE + minGres >= max_E)
				continue;
			for(int n3 = n2; n3 <= (window_size - n1 - n2)/2; n3++){
				if(Lst_OBE[n1].begin()->get_energy() + Lst_OBE[n2].begin()->get_energy() + Lst_OBE[n3].begin()->get_energy() + min_RE + minGres >= max_E)
					continue;				
				for(int n4 = n3; n4 <= (window_size - n1 - n2-n3); n4++){
					if(Lst_OBE[n1].begin()->get_energy() +Lst_OBE[n2].begin()->get_energy() + Lst_OBE[n3].begin()->get_energy() + Lst_OBE[n4].begin()->get_energy() + Gres[n1+n2+n3+n4][3]>= max_E)
                        continue;

                    // V: I added this line to calculate the residual superhelical density
                    // double a_residual = calc_alpha_res(n1+n2+n3+n4);

                    for(i = Lst_OBE[n1].begin(); i != Lst_OBE[n1].end(); ++i){
                        int p1 = i->get_pos1();
                        double e1 = i->get_energy();
                        if (e1 >= 10000) break;
                        j = Lst_OBE[n2].begin();
                        if(n1 == n2) ++j; // no repeat of the same group
                        if(e1 + j->get_energy() + Lst_OBE[n3].begin()->get_energy() + Lst_OBE[n4].begin()->get_energy() + Gres[n1+n2+n3+n4][3] >= max_E) break;
                        for(; j != Lst_OBE[n2].end(); ++j){
                            int p2 = j->get_pos1();
                            double e2 = j->get_energy();
                            if (e2 >= 10000) break;
                            k = Lst_OBE[n3].begin();
                            if(n2 == n3) ++k;
                            if(e1 + e2 + k->get_energy() + Lst_OBE[n4].begin()->get_energy() + Gres[n1+n2+n3+n4][3] >= max_E) break;
                            for(; k != Lst_OBE[n3].end(); ++k){
                                int p3 = k->get_pos1();
                                double e3 = k->get_energy();
                                if (e3 >= 10000) break;
                                l = Lst_OBE[n4].begin();
                                if(n3 == n4) ++l;
                                if(e1 + e2 + e3 + l->get_energy() + Gres[n1+n2+n3+n4][3] >= max_E) break;
                                for(; l != Lst_OBE[n4].end(); ++l){
                                    int p4 = l->get_pos1();
                                    double e4 = l->get_energy();
                                    if (e4 >= 10000) break;
                                    double e = e1 + e2 + e3 + e4 +Gres[n1+n2+n3+n4][3];
                                    if(e >= max_E) break;
                                    // checking for overlaps
                                    if(!overlap(p1, p2, n1, n2) && !overlap(p1, p3, n1, n3) && !overlap(p1, p4, n1, n4) && !overlap(p2, p3, n2, n3) && !overlap(p2, p4, n2, n4) && !overlap(p3, p4, n3, n4)){
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
                                            update_promatrix(p4, n4, e, exponent);
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
	StatesInFourRuns = (long)count;
    ZsumB += sum_4RB;
	ZsumG += sum_4RG;

    // V: Scaling factor included
    double prob = (ZsumB-exp(-alpha*alpha*K/2/RT + scaling_factor))/ZsumB;
    
	StatesInFourRuns = (long)count;
	totalFreq = StatesInOneRun + StatesInTwoRuns + StatesInThreeRuns+StatesInFourRuns;
    double runs= (sum_1RB+2*sum_2RB+3*sum_3RB+4*sum_4RB)/ZsumB;
    if (results) {
        cout << "Number of four-run states = " << count << endl;
        cout << "Total number of states = " << totalFreq << endl;
        // cout << "Total Partition Function = " << ZsumB << endl;
        // V: We need to multiply by the scaling factor to get the actual partition funciton.
        cout << "Total Partition Function = " << exp(-scaling_factor)*ZsumB << endl;
        // V: And just in case I'll print the scaling factor.
        cout << "Scaling factor " << scaling_factor << endl;
        cout << "Average number of runs = " << runs << endl;
        cout << "Transition Probability = " << prob << endl;

    }
	totalFreq = StatesInOneRun + StatesInTwoRuns + StatesInThreeRuns+StatesInFourRuns;
	flag_minE4 = 0;
	return flag_minE4;
}

#endif
