// qsidd.cpp: main() function
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

#include "SIDD_4R.h"
#include <time.h>
#include <iostream>
#include <unistd.h>
#include <sstream>
#include <string>
#include <stdlib.h>

// #include <list>
// #include <vector>


static char help[] = "\
QSIDD Help\n\
A program for computing the competition of melting, Z-DNA, and cruciform transitions in superhelical DNA.\n\n\
Required input:\n\
Input must include a sequence.\n\
Single sequences only, not databases.\n\
Sequence file or string must the last argument on the command line. \n\n\
Notes on options:\n\
-f  must be followed by a sequence file.  Otherwise type sequence string on the command line.\n\
-e  allows you to set an energy values for particular base pairs.  The energy assignment file has a very simple format. You specify two values:\n\
a coordinate, and an energy value. For example:\n\
123\t10\n\
124\t10\n\
125\t10\n\n\
Notes on Output:\n\
The output can display various runtime statistics (option -r) and model parameters (option -p).\n\
By default the output is a table showing the position, probability of melting, Z-DNA, and cruciform transitions\n\
for each nucleotide in the sequence.\n\
The nucleotide can be added as a column with the -b option.\n\
Example of first six lines of output:\n\
Position    P_melt  P_Z     P_cruciform\n\
1   2.38179e-07     2.64249e-05     0\n\
2   3.03874e-07     8.46169e-06     0\n\
3   2.83222e-07     1.19353e-05     0\n\
4   4.34397e-07     8.6928e-06      0\n\
5   4.60069e-07     4.63088e-05     0\n\
";

static char usage[] = "\
usage: qsidd [options] <sequence> \n\
options:\n\
-h        help\n\
-f        <file> sequence\n\
-c        circular DNA [default linear]\n\
-f        set if sequence in file format. STDIN assumed otherwise. This must be the last argument! \n\
-n        nearest neighbor energetics [default copolymer]\n\
-e        <file> energy assignment\n\
-i  	  set ionic strength [default 0.01]\n\
-T        set temperature [default 310.00]\n\
-s        set stress level [default -0.06]\n\
-t        set energy threshold [default 12]\n\
-X        a string containing IR positions, lengths, and energies obtained from IR-finder\n\
-b        display base in output \n\
-p        print algorithm parameters\n\
-r        print ensemble average results\n\
";


int main(int argc, char* argv[])
{
	SIDD_4R sidd;
	time_t time_1, time_2;
	int c;
	extern int optind;
	
	// default parameters
	Energetics et = Copolymeric;
	Molecule mt   = Linear;
	char *energy  = NULL;
	char *dnafile = NULL;

    stringstream dnasequence;
    stringstream cruciform_string;
    int showbase = 0;
	int verbose = 0;
    int showpar = 0;
    int showres = 0;
	double salt_conc = 0.01;
	double temperature = 310.00;
	double stress_level=0.06; 
    double threshold = 12;
	int usefile = 0;
	
    // option processing
    while ((c = getopt(argc, argv, "hcfnbsvpriTXte:")) != -1) {
			switch (c) {
			case 'c': mt = Circular;      break;
			case 'n': et = Near_Neighbor;      break;
            case 'f': usefile = 1; break;
            case 'e': energy = argv[optind++];    break;
            case 'b': showbase = 1; break;
            case 'v': verbose = 1;  break;
            case 'p': showpar = 1;  break;
            case 'r': showres = 1;  break;
			case 'T':temperature=atof(argv[optind++]);break;
			case 's':stress_level=atof(argv[optind++]);break;
            case 'i': salt_conc=atof(argv[optind++]);break;
            case 't': threshold=atof(argv[optind++]);break;
            case 'X': cruciform_string << argv[optind++]; break;	
			case 'h': cout<<help<<endl;   exit(1);
			default:                      exit(1);
		}
	}
    //start time
	time_1 = time(NULL); 
    sidd.set_showres(showres);
	sidd.set_salt_conc(salt_conc);
	sidd.set_stress_level(stress_level);
	sidd.set_temperature(temperature);
	sidd.set_MoleculeType(mt);                    
	sidd.set_EnergyType(et);
    sidd.set_threshold(threshold);
    sidd.set_cruciform_string(cruciform_string.str().c_str());
    sidd.set_showbase(showbase);

	if (argc - optind != 1) {
		cerr << usage << endl;
		exit(1);
	}
    	
	if(usefile){
        dnafile = argv[argc -1 ];
        if (showpar)
            cout << "DNA sequence file: " << dnafile << endl;
		ifstream f(dnafile);
        if (f){
            dnasequence << f.rdbuf();    
            f.close();
        }   
    }   
    else {
        dnasequence << argv[argc -1];  
    }
    int length = dnasequence.str().size();
    
	
	if(verbose)  
	{
        if (energy) {
            cout << "Energy assignment file: " << energy << endl;
        }	
        sidd.set_Flag_PEA(energy);
	}
	
    if (!sidd.initializer(dnasequence.str().c_str(),length)) {
		cerr << "sidd.initializer failed" << endl;
		exit(1);
	}
    if (sidd.get_Flag_PEA()) sidd.assign_pea(energy);
    
    if (threshold < 9)
        cout << "WARNING: threshold is too small, results may be inaccurate"<< endl;
    if (threshold > 15)
        cout << "WARNING: threshold is too high, execution time may be very long"<< endl;
    if (stress_level > 0.15 || stress_level < -0.15)
        cout << "WARNING: superhelical density is outside of physiological range"<< endl;
    if (temperature < 220 || temperature > 320)
        cout << "WARNING: temperature is outside of physiological range"<< endl;
    if (salt_conc < 0.0001)
        cout << "WARNING: salt concentration is outside of physiological range"<< endl;
    if (length < 1500)
        cout << "WARNING: sequence length is too short"<< endl;
    if (length > 10000)
        cout << "WARNING: sequence length is too long"<< endl;

    //sidd.show_seq();      //	show sequence for debugging

    // V: In the code, they use the word Open very often, but it doesn't always refer to melting, it is more like the
    //    calculation of the energy at a given state (melt, Z, cruci...). It seems that it was recycled from the
    //    strand-separation code SIDD.
    // V: gen_OpenBaseEnergy() actually calculates the energies of states with 1 run for each type of transition.
    //    I think the objective is to obtain the 2-dimensional vector/list Lst_OBE[window_size][transition].
    //    Each entry of this vector, is a list containing objects of type stat_1R with associated position and energy.

	sidd.gen_OpenBaseEnergy(); // SIDD_1R (calculate opening energies and sort w/ increasing energy for each window size 1 to 250)
	sidd.write_close(); //don't write profile yet, just find minE
	sidd.Search_Low1RE(); // this step is here to update minE (from the zero-run state) if it's found in one-run states
	sidd.write_open(); // start to store info
	if(verbose)
        cout << "writing profile...\n";
	if (showpar)
        sidd.show_parameter();
	sidd.reset_profile();  // initialize profile
	sidd.Search_Low1RE(); // searching for one-run states and now storing info 
	sidd.Search_Low2RE(); // searching states for two-run
	sidd.Search_Low3RE(); // searching states for three-run
	sidd.Search_Low4RE(); // searching states for four-run

    // V: I think for all these Search_LowNRE, global quantities such as the number of states, partition function,
    // energies Boltzmann factors and everything is calculated (the method is applied).
    // V: But the following lines actually calculate parameters for each region on the DNA
    // (e.g., the probability of transitioning each bp).

    // V: I am not sure what fill_profile() does, but I think it prepares the array with the profiles...
	sidd.fill_profile(); // computing profile
    // V: calc_profile does what it says, it calculates the profiles...
	sidd.calc_profile(); // calculate profile
    // V: This one sums probabilites, which basically gives you the number of runs in a particular transition.
    sidd.sum_open();
    time_2 = time(NULL);  //end time
    if (showpar)
        cout << "Run time = " << (time_2-time_1) << " sec" << endl;
    // V: This is where they print everything.
	sidd.show_profile(); // send output to screen or a disk file
	
	return 0; // end of program
    
}
