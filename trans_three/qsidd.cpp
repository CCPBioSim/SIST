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

static char help[] = "\
QSIDD Help\n\
A program for computing probabilities of melting, Z-DNA, or cruciform transitions in superhelical DNA.\n\
To analyze the cruciform transition, you must first run IR_finder.pl.\n\n\
Required input:\n\
Input must include a sequence.\n\
Single sequences only, not databases.\n\
Sequence file or string must the last argument on the command line.\n\
Melting transition is the default.\n\n\
Notes on options:\n\
-f  must be followed by a sequence file.  Otherwise type sequence string on the command line.\n\
-Z  set for Z-DNA transition.\n\
-C  set for cruciform transition, must also include -X followed by a string produced by IR_finder.pl\n\
-e  allows you to set an energy values for particular base pairs.  The energy assignment file has a very simple format. You specify two values:\n\
a coordinate, and an energy value. For example:\n\
123\t10\n\
124\t10\n\
125\t10\n\n\
Notes on Output:\n\
The output can display various runtime statistics (option -r) and model parameters (option -p).\n\
By default the output is a table showing the position, probability of specified transition for each nucleotide in the sequence.\n\
For the melting transition the energy of destabilization G(x) is also displayed.\n\
The nucleotide can be added as a column with the -b option.\n\
Example of first six lines of output for the melting transition:\n\
Position        P(x)    G(x)\n\
1       3.88666e-07     11.4338\n\
2       4.96282e-07     11.4515\n\
3       4.63064e-07     11.4661\n\
4       7.09681e-07     11.4445\n\
5       7.51455e-07     11.4479\n\
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
  -i  	    set ionic strength [default 0.01]\n\
  -T 	    set temperature [default 310.00]\n\
  -s        set stress level [default -0.06]\n\
  -t        set energy threshold [default 12]\n\
  -Z        Z-DNA calculation with dinucleotide energetics [default melting]\n\
  -C        Indicates the cruciform calculation will be performed [default melting]\n\
  -X        a string containing IR positions, lengths, and energies obtained from IR-finder\n\
  -b        include base in result file\n\
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

	int verbose = 0;
    int showpar = 0;
    int showres = 0;
    int showbase = 0;
	double salt_conc = 0.01;
	double temperature = 310.00;
	double stress_level=0.06; 
	double threshold = 12;
	int usefile = 0;
	// option processing
	while ((c = getopt(argc, argv, "hcfnbvpraietTZCXsm:")) != -1) {
			switch (c) {
			case 'c': mt = Circular;    break;
			case 'n': et = Near_Neighbor;   break;
            case 'b': showbase = 1; break;
            case 'f': usefile = 1; break;
			case 'Z': et = Z_DNA;   break;
            case 'C': et = Cruciform;   break;
            case 'X': cruciform_string << argv[optind++];   break;
			case 'e': energy = argv[optind++];  break;
			case 'v': verbose = 1;  break;
            case 'p': showpar = 1;  break;
            case 'r': showres = 1;  break;
			case 'T': temperature = atof(argv[optind++]);   break;
			case 's': stress_level = atof(argv[optind++]);  break;
			case 'i': salt_conc=atof(argv[optind++]);   break;
			case 't': threshold=atof(argv[optind++]);   break;
			case 'h': cout<<help<<endl;   exit(1);
            default:  exit(1);
		}
	}
	
	//start time
	time_1 = time(NULL); 
	sidd.set_salt_conc(salt_conc);
	sidd.set_stress_level(stress_level);
	sidd.set_temperature(temperature);
	sidd.set_MoleculeType(mt);                    
	sidd.set_EnergyType(et);
	sidd.set_threshold(threshold);
    sidd.set_showres(showres);
    sidd.set_cruciform_string(cruciform_string.str().c_str());
    sidd.set_showbase(showbase);
    
	if (argc - optind != 1) {
		cerr << usage << endl;
		exit(1);
	}

	if(usefile){
        dnafile = argv[argc -1 ];
        if (showpar) {
            cout << "DNA sequence file: " << dnafile << endl;
        }
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

	if(verbose)  {
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

//	sidd.show_seq(); // show sequence
	sidd.gen_OpenBaseEnergy(); // SIDD_1R (calculate opening energies and sort w/ increasing energy for each window size 1 to 250)
	sidd.write_close(); //don't write profile yet, just find minE
	sidd.Search_Low1RE(); // this step is here to update minE (from the zero-run state) if it's found in one-run states
	sidd.write_open(); // start to store info
	if(verbose)
        cout << "writing profile...\n";
    if (showpar)
        sidd.show_parameter();
    //cout << "AAAA1"<< endl;
	sidd.reset_profile();  // initialize profile
    //cout << "AAAA2"<< endl;
	sidd.Search_Low1RE(); // searching for one-run states and now storing info
    //cout << "AAAA3"<< endl;
	sidd.Search_Low2RE(); // searching states for two-run
    //cout << "AAAA4"<< endl;
	sidd.Search_Low3RE(); // searching states for three-run
    //cout << "AAAA5"<< endl;
	sidd.Search_Low4RE(); // searching states for four-run
    //cout << "AAAA6"<< endl;
	sidd.fill_profile(); // computing profile
    //cout << "AAAA7"<< endl;
	sidd.calc_profile(); // calculate profile
	time_2 = time(NULL);  //end time
    sidd.sum_open();
    time_2 = time(NULL);  //end time
    if (showpar)
        cout << "Run time = " << (time_2-time_1) << " sec" << endl;
    sidd.show_profile(); // send output to screen or a disk file	
	return 0; // end of program

}
