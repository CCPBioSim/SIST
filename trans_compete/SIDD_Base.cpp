// SIDD_Base.cpp: implementation of the SIDD_Base class.
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

#ifndef SIDD_BASE_CPP
#define SIDD_BASE_CPP

#include "SIDD_Base.h"
#include <stdlib.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SIDD_Base::SIDD_Base()
{
	plasmid_seq = 0; // null pointer
    en_cruciforms = 0;
	for (int t=0; t<=Ns; t++) {
		profile[t] = 0;
	for(int m = 0; m <= MaxInitialWindowSize; m++) {
		promatrix[m][t] = 0;
		}
	}
    
	R = 8.314/(4.2*1000.0); // gas constant
    C = 3.6; // stiffness coefficient
	A = 10.4; // bps per turn in B-DNA
	Az = -12.0; //bps per turn in Z-DNA
	tz = 0.4; //undertwist at the two junctions of Z-DNA
	alpha = 0.0; // linking difference
	K0 = 2220.0; // linking energy coefficient	
	pea = 0; // null pointer
	Flag_PEA = false;
	EnergyType = Copolymeric; // default value
	MoleculeType = Linear; // default value

	min_E = 0.0; // minimum energy
	max_E = 0.0;
	min_WS = 0;  // the window size corresponding to the minimum energy
	write_profile = false; // flagging if writing starts

	ZsumB = 0.0; // summation of Boltzman frequence
	ZsumG = 0.0; // store partition function value
	sum_1RG = sum_2RG = sum_3RG = 0.0;
	sum_1RB = sum_2RB = sum_3RB = 0.0;
    // MEEE: I ADDED THIS FOR CALCULATING NUUMBER OF RUNS ACCURATELY
    sum_1RB = sum_2RB = sum_3RB = 0.0;
    set_showbase(1);
}

SIDD_Base::~SIDD_Base()
{
    
	delete [] plasmid_seq;
	delete [] sa; 
	delete [] pea;
    for(int i = 0; i < length_seq; i++) {
        delete [] en_cruciforms[i];
    }
	for (int t=0; t<=Ns; t++) {
		delete [] profile[t]; 
		for(int i = MinWindowSize; i <= MaxWindowSize; i++) {
            delete [] promatrix[i][t];
		}
    }
}

void SIDD_Base::set_stress_level(double stress_level)
{
	if(stress_level > 0)stress_level = -stress_level;
	supdensity = stress_level; 
}

void SIDD_Base::set_threshold(double threshold)
{
	theta = threshold;
}

//energy parameters
void SIDD_Base::set_salt_conc(double salt_conc)
{
	Salt_Conc = salt_conc;
	TMAT = 354.65 + 16.6*log10(Salt_Conc); 
	TMGC = TMAT + 41.0;
}

//energy parameters
void SIDD_Base::set_temperature(double temperature)
{
	Temperature = temperature;
	RT = 1.9872*Temperature/1000.0; // constant
 	BAT = 7.2464*(1.0 - Temperature / TMAT); // coefficient (kcal)
	BGC = 9.0172*(1.0 - Temperature / TMGC); // coefficient (kcal)
}

void SIDD_Base::set_Cinitiation()
{
    Ecr = 192.5-Temperature*0.565-4*BAT-2*2.44*RT*log(4);
}


void SIDD_Base::set_EnergyType(Energetics et)
{
	EnergyType = et;
}

void SIDD_Base::set_min_e(double min_e)
{
          min_E = min_e;
          max_E = min_E + theta;
}
void SIDD_Base::set_MoleculeType(Molecule mt)
{
	MoleculeType = mt;
}

void SIDD_Base::set_MaxWindowSize()
{
	MaxWindowSize = 350;
	if (length_seq <  MaxWindowSize){
		MaxWindowSize = length_seq;
	}	
}

void SIDD_Base::set_MinWindowSize()
{
	MinWindowSize = 1;
}

// linking energy coefficient
void SIDD_Base::set_K()
{
	K = K0*RT / length_seq;
}

void SIDD_Base::set_cruciform_string(std::string cruciform_string)
{
	cr_string = cruciform_string;
} 

//nucleation energies
void SIDD_Base::set_junction()
{
  a[0]=10.84;  // Melt junction energy
  a[1]=10.0; // Z-DNA junction energy	
  a[2]=0.0; //zero for cruciform since the energy is included in IR input string 
}

//linking difference
void SIDD_Base::set_Alpha()
{
	alpha = supdensity*length_seq / A;

}

//sequence length
int SIDD_Base::get_sequence_length(){
	if(MoleculeType == Linear)
		return (length_seq - Len_Gap_Seq);
	return length_seq;
}

void SIDD_Base::set_sequence_length(int len){
	length_seq = len;
	if(MoleculeType != Circular)
        length_seq+=Len_Gap_Seq;
}

void SIDD_Base::set_showres(int showres)
{
	results = showres;
}

//read sequence
double SIDD_Base::prepare_sequence(std::string sequence) {
    plasmid_seq = new int[length_seq+1];    
    for(int j = 0; j < length_seq; j++){
        plasmid_seq[j] = 0;
    }
    int read = 1;
    int i = 0;	
    string::iterator iter;
    for(iter=sequence.begin();iter !=sequence.end();iter++)
    {    
        char ch = *iter;	
        if(ch == '>') {
            read = 0;
        }
        if(!read) {
            if(ch == '\n')  
                read = 1;
            continue;
        }
        if(read) {
            ch = toupper(ch);
            switch(ch)
            {
                case 'A':
                    plasmid_seq[i++] = 0;
                    break;
                case 'C':
                    plasmid_seq[i++] = 1;
                    break;
                case 'G':
                    plasmid_seq[i++] = 2;
                    break;
                case 'T':
                    plasmid_seq[i++] = 3;
                    break;
                //turn Ns into Gs
                case 'N':
                    plasmid_seq[i++] = 2;
                    break;
                default: 
                    break;
            }
        }
    }
    length_seq = i;
    if(!length_seq) {
        cerr << "sequence file not available.\n";
        return false;
    }
	
  	if(MoleculeType != Circular)
		length_seq+=Len_Gap_Seq;
	return true;
}

//read IR input string
void SIDD_Base::prepare_cruciforms() {
    string str1 = cr_string;
    size_t found1 = 0;
    int IR[2] = {0, 0};
    double energy = 0.0;
    int pos1 = 0;
    int pos2 = 0;
    string str2;
    string str3;
    
    en_cruciforms = new double* [length_seq+1];    
    for(int i = 0; i < length_seq; i++){
        en_cruciforms[i] = new double[MaxWindowSize+1];
    }
    
    for(int i = 0; i < length_seq; i++){
        for(int j = 1; j <= MaxWindowSize; j++){
            en_cruciforms[i][j] = 10000;
        }
    }
    
    while (found1!=string::npos) {
        found1=str1.find("|",pos1);
        if(found1!=string::npos) {
            str2 = str1.substr(pos1,found1-pos1); 
            pos2 = 0;
            size_t found2 = 0;
            int i = 0;
            while (found2!=string::npos) {
                found2=str2.find(",",pos2);
                if(found2!=string::npos) {
                    str3 = str2.substr(pos2,found2-pos2);
                    pos2 = int(found2)+1;
                    if (i == 2)
                        energy = atof(str3.c_str());
                    else if (i < 2)
                        IR[i] = atoi(str3.c_str());
                    i++;
                }
            }
            if (i == 3 && IR[0] > 0 && IR[0] <= length_seq &&
                IR[1] > 0 && IR[1] <= MaxWindowSize)
                en_cruciforms[IR[0]-1][IR[1]] = energy;
            pos1 = int(found1)+1;
        }
    }
 }

 //TODO: This is where the energetics are stored. But what does it actually mean each of the temrs...?
//nearest neighbor energetics
void SIDD_Base::set_Delta_G()
{
	double Delta_H[4][4] = {
		{32.3, 35.8, 36.0, 31.2},
		{35.8, 35.1, 39.8, 36.0},
		{36.0, 39.8, 35.1, 35.8},
		{31.2, 36.0, 35.8, 32.3}
	};
    
	double Delta_S[4][4] ={ // J per mol per Kelvin degree
		{95.6, 100.4, 102.4, 95.3},
		{100.4, 93.9, 102.4, 102.4},
		{102.4, 102.4, 93.9, 100.4},
		{95.3, 102.4, 100.4, 95.6}
	};
    
    
	// unit conversion: KJ to Kcal
	// 1 kcal = 4.184 KJ
	int i, j;
	for(i = 0; i < 4; i++){
		for(j = 0; j < 4; j++){
			Delta_H[i][j] /= 4.184;
			Delta_S[i][j] /= 4184.0;
		}
	}

	for(i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			Delta_S[i][j] = 1.0/(16.6*log10(Salt_Conc/0.1)/Delta_H[i][j] + 1.0 / Delta_S[i][j]);
		}
	}
    
	for(i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
			Delta_G[i][j] = Delta_H[i][j] - Temperature*Delta_S[i][j]; // kcal per mole
		}
	}
}

//superhelical energy
double SIDD_Base::calc_Gres(int n_m, int n_z, int n_c, int nr)
{
    // V: This has the form of 2013 paper competitive transitions DIna and Craig.
    double  t1 = (alpha + n_m/A+ n_z/A - n_z/Az + 2*tz*nr + n_c/A) * (alpha + n_m/A+ n_z/A - n_z/Az + 2*tz*nr + n_c/A);	
	double t2 = 2*PI*PI*C*K*t1 /(4*PI*PI*C + K*n_m);
	return t2;
}

void SIDD_Base::set_E_limit()
{
	min_E = 0.5*K*alpha*alpha;
	max_E = min_E + theta;
}

//minimum residual energy
void SIDD_Base::set_minGres()
{    
minGres = calc_Gres(1,0,0,0);
	for (int r = 0; r <= 4; r++) {
        for(int w1 = 0; w1 <= MaxWindowSize; w1++){
            for(int w2 = 0; w2 <= MaxWindowSize; w2++){
                for(int w3 = 0; w3 <= MaxWindowSize; w3++){
                    if(calc_Gres(w1,w2,w3,r) < minGres)
                        minGres=calc_Gres(w1,w2,w3,r);
                }
            }
        }
    }
}

int SIDD_Base::delta_fnc(int x, int y) {
    return (x == y);
}

// computing transition energy: Z-DNA
void SIDD_Base::calc_Z(int pos)
{
    
    double energy_Z_AS[4][4] = {
		{3.9,4.6,3.4,5.9},
		{1.3,2.4,0.7,3.4},
		{3.4,4.0,2.4,4.6},
		{2.5,3.4,1.3,3.9}
	}; 
	
	double energy_Z_SA[4][4] = {
		{3.9,1.3,3.4,2.5},
		{4.6,2.4,4.0,3.4},
		{3.4,0.7,2.4,1.3},
		{5.9,3.4,4.6,3.9}
	}; 
	
	double energy_Z_ZZ[4][4] = {
		{7.4,4.5,6.3,5.6},
		{4.5,4.0,4.0,6.3},
		{6.3,4.0,4.0,4.5},
		{5.6,6.3,4.5,7.4}
	}; 
 
	// Base code before pos
	int p1 = plasmid_seq[pos-1 < 0 ? length_seq - 1 : pos - 1]; 
	// Base code after pos
	int p2 = plasmid_seq[pos+1 == length_seq ? 0 : pos + 1];
	// Base code of pos 
	int p = plasmid_seq[pos]; 
	// Config code (anti or syn) before pos
	int c1 = sa[pos-1 < 0 ? length_seq - 1 : pos - 1];
	// Config code (anti or syn) after pos
	int c2 = sa[pos+1 == length_seq ? 0 : pos + 1];
	// Config code (anti or syn) of pos
	int c = sa[pos];
	e1=0;
	e2=0;
	zz=0;
	if (c != c2) 
	{
		if (c == 'a')	{
			e1 = energy_Z_AS[p][p2];
		}
		else {
			e1 = energy_Z_SA[p][p2];
		}
	}
	
	if (c != c1) 
	{
		if (c == 'a')	{
			e2 = energy_Z_AS[p][p1];
		}
		else {
			e2 = energy_Z_SA[p][p1];
		}
	}
	if (c == c1) {
		zz =energy_Z_ZZ[p][p1];
		int p0 = sa[pos-2 < 0 ? length_seq - 2 : pos - 2];
		if (energy_Z_ZZ[p][p2] >= energy_Z_ZZ[p0][p1]) {
			if (c=='a')
				e2 = energy_Z_SA[p1][p];
			else 
				e2 = energy_Z_AS[p1][p];
		}
		else {
			if (c1=='a')
				e2 = energy_Z_AS[p1][p];
			else 
				e2 = energy_Z_SA[p1][p];
		}
		
	}	
	if (c == c2) {
		int p3 = plasmid_seq[pos+2 == length_seq ? 1 : pos + 2];
		if (energy_Z_ZZ[p2][p3] >= energy_Z_ZZ[p][p1]) {
			zz = energy_Z_ZZ[p][p1]; 
			if (c2=='a')
				e1 = energy_Z_SA[p][p2];
			else 
				e1 = energy_Z_AS[p][p2];
		}
		else {
			zz = energy_Z_ZZ[p2][p3]; 
			if (c=='a')
				e1 = energy_Z_AS[p][p2];
			else 
				e1 = energy_Z_SA[p][p2];
		}
	}		
}


// start position - startp
// window size - n
// return - sum of window energy
double SIDD_Base::sum_WindowZ(int startp, int n)
{
	double s = 0.0;
	if (n%2==1 || n==2 || n==4 || n==6) {
		s=10000;
		return s;
	}
	else { 
		if(startp + n - 1 < length_seq){
			for(int i = startp; i < startp + n; i++) {
				calc_Z(i);
				if (i==startp) 
					s+=e1/2;
				if ((i-startp)%2==0 && i!=startp)
					s+=e1/2+zz;
				if ((i-startp)%2==1)
					s+=e2/2;
			}
			return s;
		}
		else{
			for(int i = startp; i < length_seq; i++) {
				calc_Z(i); 
				if (i==startp) 
					s+=e1/2;
				if ((i-startp)%2==0 && i!=startp)
					s+=e1/2+zz;
				if ((i-startp)%2==1)
					s+=e2/2;
			}
			for(int j = 0; j < startp + n - length_seq; j++) {
				calc_Z(j);
				if (j==startp) 
					s+=e1/2;
				if ((length_seq-startp+j)%2==0 && j!=startp)
					s+=e1/2+zz;
				if ((length_seq-startp+j)%2==1)
					s+=e2/2;
			}
			return s;
		}
	}
}


// computing transition energy: neighbor interation
double SIDD_Base::calc_NI(int pos)
{
	if(pos < length_seq && pos >= 0){

		if(Flag_PEA && pea[pos] != 0.0)
			return pea[pos]; // assign specified energy
		int p1 = plasmid_seq[pos-1 < 0 ? length_seq - 1 : pos - 1]; 
		int p2 = plasmid_seq[pos+1 == length_seq ? 0 : pos + 1]; 
		int p = plasmid_seq[pos]; 
		return (Delta_G[p1][p] + Delta_G[p][p2]) / 2.0;
	}
	else{
		cerr << "position: out of range.\n";
		return 0;
	}
}

// start position - startp
// window size - n
// return - sum of window energy
double SIDD_Base::sum_WindowNI(int startp, int n)
{
	double s = 0.0;
	if(startp + n - 1 < length_seq){
		for(int i = startp; i < startp + n; i++)
			s += calc_NI(i);
		return s;
	}
	else{
		for(int i = startp; i < length_seq; i++)
			s += calc_NI(i);
		for(int j = 0; j < startp + n - length_seq; j++)
			s += calc_NI(j);
	}
	return s;
}

double SIDD_Base::sum_WindowNIbyBases(int startp, int n)
{
	int at, gc;
	at = gc = 0;
    double  s = 0.0;
	double tempE = 0.0;
    count_AT_GC(startp, n, at, gc, tempE);
    s = BAT*at + BGC*gc + tempE;
    return s;
}

double SIDD_Base::sum_cruciform(int startp, int n)
{
    return en_cruciforms[startp][n];
}

// counting A/T bases
void SIDD_Base::count_AT_GC(int startp, int n, int& c_AT, int& c_GC, double& tempE)
{
	c_AT = c_GC = 0;
	tempE = 0.0;
	if(startp + n - 1 < length_seq){
		for(int i = startp; i < startp + n; i++){
			if(Flag_PEA && pea[i] != 0.0)
				tempE += pea[i];
			else{
				if(plasmid_seq[i] == 0 || plasmid_seq[i] == 3)
					c_AT++;
				else
					c_GC++;
			}
		}
	}

	else{
		for(int i = startp; i < length_seq; i++){
			if(Flag_PEA && pea[i] != 0.0)
				tempE += pea[i];
			else{

				if(plasmid_seq[i] == 0 || plasmid_seq[i] == 3)
					c_AT++;
				else
					c_GC++;
			}
		}

		for(int j = 0; j < startp + n - length_seq; j++){
			if(Flag_PEA && pea[j] != 0.0)
				tempE += pea[j];
			else{
				if(plasmid_seq[j] == 0 || plasmid_seq[j] == 3)
					c_AT++;
				else
					c_GC++;
			}
		}
	}
}


// energy of a transition for a sequence segment
// starting at starp of length n for transition defined by type
double SIDD_Base::calc_OPenBasesEnergy(int startp, int n, int type)
{
	if (type == 0) {
        switch(EnergyType){
            case Near_Neighbor:
                return sum_WindowNI(startp, n);
            case Copolymeric:
                return sum_WindowNIbyBases(startp, n);
            default:
                return sum_WindowNIbyBases(startp, n);
        }
    }
    else if (type == 1) 
		return sum_WindowZ(startp,n);
    else
        return sum_cruciform(startp,n);        
}


// reading sequence and initializing parameters
bool SIDD_Base::initializer(std::string sequence,int len)
{
    set_sequence_length(len);
    if(!prepare_sequence(sequence)) {
		//sequence either 0 or too long
		return 0;
	}
    	set_MaxWindowSize();
	set_MinWindowSize();
	set_K();
	set_Alpha();
	set_junction();
	set_Delta_G();
	set_E_limit();
	prepare_cruciforms();
    
	sa = new char[length_seq+1];
	pea = new double[length_seq+1];
	for(int i = 0; i < length_seq+1; i++)
		pea[i] = 0.0;

	for (int t=0; t<=Ns; t++) {
		profile[t] = new G_x[length_seq+1];
	for(int m = MinWindowSize; m <= MaxWindowSize; m++) {
		promatrix[m][t] = new G_x[length_seq+1];
		}
	}
	
	if(MoleculeType == Linear){
		for(int k = length_seq - Len_Gap_Seq; k < length_seq; k++)
			plasmid_seq[k] = 2;  //add a tail of G's
	}
    
	set_minGres();
	reset_profile();

    for(int i = 0; i < length_seq; i++){			
        switch(plasmid_seq[i]){
				case 0: 
					sa[i]='s';
					break;
				case 1: 
					sa[i]='a';
					break;
				case 2: 
					sa[i]='s';
					break;
				case 3: 
					sa[i]='a';
					break;
			}
		}
    for(int i = 0; i < length_seq; i++){
			// Config code (anti or syn) before pos
			int c1 = sa[i-1 < 0 ? length_seq - 1 : i - 1];
			// Config code (anti or syn) after pos
			int c2 = sa[i+1 == length_seq ? 0 : i + 1];
			// Config code (anti or syn) of pos
			int c = sa[i];
			if (c==c1 && c==c2) {
				if (c=='s')
					sa[i]='a';
				else
					sa[i]='s';	
			}
    }
		
    for(int i = 0; i < length_seq; i++){
			// Config code (anti or syn) after pos
			int c2 = sa[i+1 == length_seq ? 0 : i + 1];
			// Config code (anti or syn) of pos
			int c = sa[i];	
			int c3 = sa[i+2 == length_seq ? 1 : i + 2];
			int c4 = sa[i+3 == length_seq ? 2 : i + 3];
			int c5 = sa[i+4 == length_seq ? 3 : i + 4];
			if (c==c2 && c3==c4 && c!=c3) {
				if (c5=='a') {
					sa[i]='a';
					sa[i+1]='s';
					sa[i+2]='a';
					sa[i+3]='s';
				}
				else {
					sa[i]='s';
					sa[i+1]='a';
					sa[i+2]='s';
					sa[i+3]='a';
				}
			}
		}
    
 return true;
}

// assign energy position-wide
void SIDD_Base::assign_pea(char* efile)
{
	ifstream ins(efile);
    
	if(ins.is_open()){
		while(!ins.eof()){
			double e0;
			int p;
			ins >> p >> e0;
			if(p >=0 && p < length_seq && e0 != 0.0){
				pea[p] = e0;
				cout << "assign energy at " << p << ": " << pea[p] << endl;
			}
		}
	}
	else{
		cerr << "fail to open the file: " << efile << endl;
		Flag_PEA = false;
	}
	ins.close();
}

void SIDD_Base::reset_promatrix()
{
	ZsumB = 0.0;
	ZsumG = 0.0;
	totalFreq = 0;
	for (int t=0; t<=Ns; t++) {
	for(int k = MinWindowSize; k <= MaxWindowSize; k++){
		for(int i = 0; i < length_seq; i++){
			promatrix[k][t][i].reset();
			}
		}
	}
}

// startp - start position
// n - window size
// x - free energy
bool SIDD_Base::update_promatrix(int startp, int n, int struc, double x, double bzfactor)
{
	if(startp < 0 || startp >= length_seq || n > MaxWindowSize)
		return false;
	
	promatrix[n][struc][startp].add(x, bzfactor, RT);
    return true;	
}

void SIDD_Base::fill_profile()
{
	for(int t=0; t<=Ns; t++) {
	for(int n = MinWindowSize; n <= MaxWindowSize; n++){
		for(int startp = 0; startp < length_seq; startp++){
			double lastG = promatrix[n][t][startp].get_sum_xG();
			double lastB = promatrix[n][t][startp].get_sum_xB();
			if(startp + n - 1 < length_seq){
				for(int i = startp; i < startp + n; i++){
					profile[t][i].add(-1.0, 0.0, -1.0, lastG, lastB); 
				}
			}
			else{
				for(int i = startp; i < length_seq; i++)
					profile[t][i].add(-1.0, 0.0, -1.0, lastG, lastB);  
				for(int k = 0; k < startp + n - length_seq; k++)
					profile[t][k].add(-1.0, 0.0, -1.0 , lastG, lastB);  
				}				
			}
		}
	}
}

void SIDD_Base::reset_profile()
{
	ZsumB = 0.0;
	ZsumG = 0.0;
	totalFreq = 0;
	for(int t=0; t<=Ns; t++) {
	for(int i = 0; i < length_seq; i++){
		profile[t][i].reset();
		}
	}

	reset_promatrix();
}


// startp - start position
// n - window size
// x - free energy
void SIDD_Base::calc_profile()
{
    // V: avg_Gs My guess is that it is the average energy.
	double ave_Gs = 0.0;
	if(ZsumB == 0.0)
		ave_Gs = 0.0;
	else
		ave_Gs = ZsumG / ZsumB;
    // V: Goes through all the bases for each transition, and calculates the free energy required for transition
    //    and the probability of transitioning.
    //    profile is of class G_x, so check those functions (calc for example).
	for(int t=0; t<=Ns; t++) {
        for(int i = 0; i < length_seq; i++){
            profile[t][i].calc(ave_Gs, ZsumB);
        }
	}
    // V: The following lines is just a search of the maximum Gx for each transition, to substitute the Infinite.
	maxGx = profile[1][0].get_ave_Gx();
	for(int t=0; t<=Ns; t++) {
        for(int j = 1; j < length_seq; j++){
            if(maxGx < profile[t][j].get_ave_Gx()){
                maxGx = profile[t][j].get_ave_Gx(); // find max Gx
            }
        }
	}
    // V: This is for the substitution.
	for(int t=0; t<=Ns; t++) {
        for(int k = 0; k < length_seq; k++){
            if(profile[t][k].get_ave_Gx() == INFINITE_Gx)
                profile[t][k].set_ave_Gx(maxGx); // replace INF Gx with maxGx
		}
	}
}

//average number of bp in each transition
void SIDD_Base::sum_open()
{
	for(int t=0; t<=Ns; t++) {
		sumn[t]=0.0;
		for(int i = 0; i < length_seq; i++) {
			sumn[t]+=profile[t][i].get_px();	
		}
	}
    if (results) {
        cout << "N_Melting = " << sumn[0] << endl;
        cout << "N_Z = " << sumn[1] << endl;
        cout << "N_Cruciform = " << sumn[2] << endl;
    }
}

void SIDD_Base::get_column_header() {
if(get_showbase()) {
    cout << "Position" << "\tBase" << "\tP_melt" << "\tP_Z"<<"\tP_cruciform"<<endl;
    // V: I added this for printing energies...
    //cout << "Position" << "\tBase" << "\tP_melt" << "\tP_Z"<<"\tP_cruciform" << "\tG_melt" << "\tG_Z"<<"\tG_cruciform"<<endl;
    }
    else {
        cout << "Position" << "\tP_melt" << "\tP_Z"<<"\tP_cruciform"<<endl;
    }
}

// TODO: V - I think here we can print th energies.
//printing output
void SIDD_Base::show_profile()
{
	if(MoleculeType == Linear)
	{
		length_seq = length_seq -Len_Gap_Seq;
	}
    get_column_header();
    for(int i = 0; i < length_seq; i++) {
        if(get_showbase()) {
            cout << i + 1 << "\t"<<decode_base(plasmid_seq[i])<<"\t"<< profile[0][i].get_px()<<"\t"<< profile[1][i].get_px()<<"\t"<<profile[2][i].get_px()<<endl;
            // V: I added this to print energies
            //cout << i + 1 << "\t"<<decode_base(plasmid_seq[i])<<"\t"<< profile[0][i].get_px()<<"\t"<< profile[1][i].get_px()<<"\t"<<profile[2][i].get_px() << "\t" << profile[0][i].get_ave_Gx()<<"\t"<< profile[1][i].get_ave_Gx()<<"\t"<<profile[2][i].get_ave_Gx() << endl;
        }
        else {
            cout << i + 1 <<"\t"<< profile[0][i].get_px()<<"\t"<<profile[1][i].get_px()<<"\t"<<profile[2][i].get_px()<<endl;
        }
    }
}

// printing free energy
void SIDD_Base::show_deltaG()
{
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++)
			cout << get_Delta_G(i,j)<< "\t";
		cout << endl;
	}
}

// printing encoded sequence
void SIDD_Base::show_seq()
{
	int k = 0;
	cout << "Total length of sequence is " << length_seq << endl;
	cout << 1 << "\t";
	for(int i = 0; i < length_seq; i++){

		switch(plasmid_seq[i]){
		case 0: 
			cout << 'A';
			break;
		case 1: 
			cout << 'C';
			break;
		case 2: 
			cout << 'G';
			break;
		case 3: 
			cout << 'T';
			break;
		default: 
			cout << 'N';
			break;
		}
		if((i+1)%10 == 0){
			cout << "  ";
			k++;
		}
		if(k%5 == 0 && (i+1)%10 == 0){
			cout << endl << i + 2 << "\t";
			k = 0;
		}
	}
	cout << "\ntotal length: " << length_seq << endl;
}

//print parameters
void SIDD_Base::show_parameter()
{
    cout << "Sequence Length = " << get_sequence_length() << endl;
    if(EnergyType == Copolymeric)
		cout << "Copolymeric melting energetics\n";
	else if(EnergyType == Near_Neighbor)
		cout << "Nearest neighbor melting energetics\n";
	if(MoleculeType == Circular)
		cout << "Molecule type: circular DNA\n";
	else
		cout << "Molecule type: linear DNA\n";
    cout << "Stress Level = " << supdensity << endl;
	cout << "linking difference = " << alpha << " turns" << endl;
	cout << "Threshold = " << theta << " kcal/mol" << endl;
	cout << "Temperature = " << Temperature << " K" << endl;
	cout << "Salt Concentration = " << Salt_Conc << " M" << endl;
	cout << "K = " << C << " kcal-bp/rad^2" << endl;
    set_Cinitiation();
    cout << "Melting nucleation energy = " << a[0] << " kcal/mol" << endl;
    cout << "Z-DNA nucleation energy = " << a[1] << " kcal/mol" << endl;
    cout << "Cruciform nucleation energy = " << Ecr << " kcal/mol" << endl;
    cout << "BAT = " << BAT << " kcal/mol; BGC = " << BGC << " kcal/mol" << endl;
	cout << "MinE = " << min_E << " kcal/mol; Max_E = " << max_E << " kcal/mol" << endl;
}

//decode sequence from numerical back to letter
char SIDD_Base::decode_base(int base)
{
    switch(base){
        case 0:  
            return 'A';
        case 1: 
            return 'C';
        case 2: 
            return 'G';
        case 3:
            return 'T';
        default:
            return 'N';
	}
}        


#endif
