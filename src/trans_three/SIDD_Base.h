// SIDD_Base.h: interface for the SIDD_Base class.
// This class is the very base one defining all the parameters and all 
// the common data and function members
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

#ifndef SIDD_BASE_H_
#define SIDD_BASE_H_
#include <list>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include "G_x.h"

const double PI = 3.1415926535897932;
const int MaxInitialWindowSize = 250;
const double GMAX_ADJ = 12;
const int Len_Gap_Seq = 50;
// V: I Added this
const double scaling_factor = 300.0;


enum Energetics{Copolymeric, Near_Neighbor, Z_DNA, Cruciform};
enum Molecule{Circular, Linear};

class SIDD_Base
{
protected:
    int showbase; //include base pairs in output
	int* plasmid_seq; // hold encoded sequence
    double** en_cruciforms; //array of cruciform pos,length,energies
	int length_seq;
	double Delta_G[4][4]; // free energy
	double Temperature;
	double Salt_Conc; // salt concentration
	double R; // constant
	double C; // torsional stiffness
	double a; // initial energy
	double A; // bases per turn
	double Az;
	double tz;
	double stresslevel;
	double supdensity;
	double alpha; // linking difference
	double theta; // threshold
	double K0;
	double K;
	double RT;
	double TMAT;
	double TMGC;
	double BAT;
	double BGC; 
	double Ecr;
	double Zmin;
	int MinWindowSize;
	int MaxWindowSize; 
	int results; 
	std::string cr_string;
	double min_E;
	double max_E;
	int min_WS;	
	char* sa;
	double e1;
	double e2;
	double zz;
	Energetics EnergyType;
	Molecule MoleculeType;
	// profile information
	G_x* profile; // hold profile
	G_x* promatrix[MaxInitialWindowSize+1];
	double minGres;
	double Gres[MaxInitialWindowSize+1][4];
	double* pea; // position-wide energy assignment
	bool Flag_PEA;
	long totalFreq;
	bool write_profile;
	double ZsumG; // total free energy
	double ZsumB; // value by partition function Z
	double sum_1RG; // summation of all states within one run
	double sum_1RB; // summation of all Boltzman factors within one run
	double sum_2RG; // summation of all states within two runs
	double sum_2RB; // summation of all Boltzman factors within two runs
	double sum_3RG;
	double sum_3RB; // summation of all Boltzman factors within three runs
	double sum_4RG;
	double sum_4RB; // summation of all Boltzman factors within three runs
	double Prob_1R;
	double Prob_2R;
	double Prob_3R;
	double maxGx;
	double minGx;

private:
	void set_K();
	void set_Delta_G();
	void set_Alpha();
	void set_junction();
	void set_E_limit();
	void set_minGres();
	void set_MinWindowSize();
	void set_MaxWindowSize();

public:
	SIDD_Base();
	virtual ~SIDD_Base();
	bool initializer(std::string filename,int len);
	void set_EnergyType(Energetics);
	void set_temperature(double temperature);
    void set_showres(int showres);
	void set_MoleculeType(Molecule);
    void set_Cinitiation();
	void show_seq();
	void show_deltaG();
	void show_parameter();
	double get_Delta_G(int i, int j){return (i< 4 && j < 4)?Delta_G[i][j]:0.0;};
	double calc_Gres(int,int);
	double calc_NI(int);
	void calc_Z(int);
	double sum_WindowNI(int startp, int n);
	double sum_WindowZ(int startp, int n);
	void count_AT_GC(int startp, int n, int& c_AT, int& c_GC, double& tempE);
	double sum_WindowNIbyBases(int startp, int n);
    double sum_cruciform(int startp, int n);
	double calc_OPenBasesEnergy(int startp, int n);
	bool update_promatrix(int startp, int n, double x, double bzfactor);
	void reset_promatrix();
	void fill_profile();
	void reset_profile();
	void get_column_header();
	void show_profile();
	void calc_profile();
    void sum_open();
	void show_avg_runs();
	void write_close(){write_profile = false;}
	void write_open(){write_profile = true;}
	void set_Flag_PEA(bool f){Flag_PEA = f;}
	bool get_Flag_PEA(){return Flag_PEA;}
	void assign_pea(char*);
	double prepare_sequence(std::string dnasequence);
    void set_cruciform_string(std::string cruciform_string);
    void prepare_cruciforms();

	void set_threshold(double threshold);
	void set_salt_conc(double salt_conc);
    void set_min_e(double min_e);
	void set_stress_level(double stress_level);
	double get_zsumb(){return ZsumB;}	
	int get_sequence_length();
	void set_sequence_length(int len);
	char decode_base(int base);
    void set_showbase(int show) {showbase = show;}
    int get_showbase() {return showbase;}
};

#endif 
