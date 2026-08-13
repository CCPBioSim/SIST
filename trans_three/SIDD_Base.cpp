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
	plasmid_seq = 0; 
	en_cruciforms = 0;
	profile = 0;
	for(int m = 0; m <= MaxInitialWindowSize; m++)
		promatrix[m] = 0;
	R = 8.314/(4.2*1000.0); // gas constant
	C = 3.6; // stiffness constant 
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
	ZsumB = 0.0; // summation of Boltzman frequencies
	ZsumG = 0.0; // store partition function value
	sum_1RG = sum_2RG = sum_3RG = 0.0;
	sum_1RB = sum_2RB = sum_3RB = 0.0;
    set_showbase(1);

}

SIDD_Base::~SIDD_Base()
{
	delete [] plasmid_seq;
	delete [] sa; 
	delete [] profile; 
	delete [] pea;
    for(int i = 0; i < length_seq; i++) {
        delete [] en_cruciforms[i];
    }
	for(int i = MinWindowSize; i <= MaxWindowSize; i++)
		delete [] promatrix[i];
}

void SIDD_Base::set_stress_level(double stress_level)
{
	if(stress_level > 0) stress_level = -stress_level;
	supdensity = stress_level; 
}

void SIDD_Base::set_threshold(double threshold)
{
	theta = threshold;
}

void SIDD_Base::set_salt_conc(double salt_conc)
{
	Salt_Conc = salt_conc;
	TMAT = 354.65 + 16.6*log10(Salt_Conc); 
	TMGC = TMAT + 41.0;
}

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
	MaxWindowSize = 250;
	if (length_seq < MaxWindowSize){ 
		MaxWindowSize = length_seq;
	}
}
void SIDD_Base::set_MinWindowSize()
{
	if (EnergyType == Z_DNA) 
		MinWindowSize = 8;
	else
		MinWindowSize = 1;	
}

void SIDD_Base::set_K()
{
	K = K0*RT / length_seq; // linking energy coefficient
}

// junction energies
void SIDD_Base::set_junction()
{
	if (EnergyType==Z_DNA) 
		a=10.0;
	else if (EnergyType==Cruciform) 
		a=0; //in string IR
	else 
        a = 10.84;
}

void SIDD_Base::set_Alpha()
{
	alpha = supdensity*length_seq/A;
}

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

void SIDD_Base::set_cruciform_string(std::string cruciform_string)
{
	cr_string = cruciform_string;
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

// reading IR string
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
        for(int j = 0; j <= MaxWindowSize; j++){
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
    profile = new G_x[length_seq+1];	
	for(int m = MinWindowSize; m <= MaxWindowSize; m++)
		promatrix[m] = new G_x[length_seq+1];
    if(MoleculeType == Linear){
		for(int k = length_seq - Len_Gap_Seq; k < length_seq; k++)
				plasmid_seq[k] = 2;  // add 50 Gs to tail when linear
	}	
	set_minGres();
	reset_profile();
	if (EnergyType == Z_DNA)
	{
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
			int c1 = sa[i-1 < 0 ? length_seq - 1 : i - 1];
			int c2 = sa[i+1 == length_seq ? 0 : i + 1];
			int c = sa[i];
			if (c==c1 && c==c2) {
				if (c=='s')
					sa[i]='a';
				else
					sa[i]='s';	
			}
		}		
		for(int i = 0; i < length_seq; i++){
			int c2 = sa[i+1 == length_seq ? 0 : i + 1];
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

//superhelical energy for each transition
double SIDD_Base::calc_Gres(int n, int nr)
{
    if (EnergyType == Z_DNA) {
        double t2 = 0.5*K*(alpha+n/A-n/Az+2*tz*nr)*(alpha+n/A-n/Az+2*tz*nr);
        return t2;
    }
    else if (EnergyType == Cruciform) {
        double  t2 = 0.5*K*(alpha + n/A) * (alpha + n/A);
        return t2;
    }
    else{
        double t1 = (alpha + n / A) * (alpha + n / A) ;
        double t2 = 2*PI*PI*C*K*t1;
        return t2 / (4*PI*PI*C + K*n);
    }
}

// V: I added this function to calculate the residual superhelical density
//    n =  number of denatured bases
//double SIDD_Base::calc_alpha_res(int n)
//{/

//    if (EnergyType == Z_DNA) {
 //       return 0.0;
 //   }
 //   else if (EnergyType == Cruciform) {
 //       return 0.0;
//    }
//    else{
 //       double t1 = 4.0*PI*PI*C*(A*alpha + n) ;
  //      double t2 = A*(4.0*PI*PI*C + n*K) ;
   //     return t1/t2 ;
 //   }
//}


void SIDD_Base::set_E_limit()
{
	min_E = 0.5*K*alpha*alpha;
	max_E = min_E + theta;
}

// set minimum residual energy
void SIDD_Base::set_minGres()
{
	minGres = Gres[MinWindowSize][0] = calc_Gres(MinWindowSize,1);
	for (int r = 0; r < 4; r++) {
        for(int w = MinWindowSize+1; w <= MaxWindowSize; w++){
            Gres[w][r] = calc_Gres(w,r+1);
            if (r == 0) {
                if(Gres[w][r] < minGres)
                    minGres = Gres[w][r];
            }
        }
    }
}

// computing transition energy: Z-DNA
void SIDD_Base::calc_Z(int pos)
{
    //Z-DNA energetics
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
	if (n%2==1) {
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

double SIDD_Base::sum_cruciform(int startp, int n)
{
    return en_cruciforms[startp][n];
}

double SIDD_Base::sum_WindowNIbyBases(int startp, int n)
{
	int at, gc;
	at = gc = 0;
	double tempE = 0.0;
	count_AT_GC(startp, n, at, gc, tempE);
	return (BAT*at + BGC*gc + tempE);

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


double SIDD_Base::calc_OPenBasesEnergy(int startp, int n)
{
	// base pair energy type
	switch(EnergyType){
	case Z_DNA:
		return sum_WindowZ(startp,n);
    case Cruciform:
        return sum_cruciform(startp,n);        
	case Near_Neighbor:
		return sum_WindowNI(startp, n);
	case Copolymeric:
		return sum_WindowNIbyBases(startp, n);
	default:
		return sum_WindowNIbyBases(startp, n);
	}

}


void SIDD_Base::reset_promatrix()
{
	ZsumB = 0.0;
	ZsumG = 0.0;
	totalFreq = 0;
	for(int k = MinWindowSize; k <= MaxWindowSize; k++){
		for(int i = 0; i < length_seq; i++){
			promatrix[k][i].reset();
		}
	}
}

// startp - start position
// n - window size
// x - free energy
bool SIDD_Base::update_promatrix(int startp, int n, double x, double bzfactor)
{
	if(startp < 0 || startp >= length_seq || n > MaxWindowSize)
		return false;
	
	promatrix[n][startp].add(x, bzfactor, RT);
//	return promatrix[n][startp].get_exponent();
 return true;	
}

void SIDD_Base::fill_profile()
{
	for(int n = MinWindowSize; n <= MaxWindowSize; n++){
		for(int startp = 0; startp < length_seq; startp++){
			double lastG = promatrix[n][startp].get_sum_xG();
			double lastB = promatrix[n][startp].get_sum_xB();
			if(startp + n - 1 < length_seq){
				for(int i = startp; i < startp + n; i++){
                    profile[i].add(-1.0, 0.0, -1.0, lastG, lastB);
				}
			}
			else{
				for(int i = startp; i < length_seq; i++)
					profile[i].add(-1.0, 0.0, -1.0, lastG, lastB);
				for(int k = 0; k < startp + n - length_seq; k++)
                    profile[k].add(-1.0, 0.0, -1.0 , lastG, lastB);
			}
		}
	}
}


void SIDD_Base::reset_profile()
{
	ZsumB = 0.0;
	ZsumG = 0.0;
	totalFreq = 0;
	for(int i = 0; i < length_seq; i++){
		profile[i].reset();
	}

	reset_promatrix();
}

void SIDD_Base::calc_profile()
{
	double ave_Gs = 0.0;
	if(ZsumB == 0.0) {
        ave_Gs = 0.0;
    }
	else{
		ave_Gs = ZsumG / ZsumB;
    }

	for(int i = 0; i < length_seq; i++){
		profile[i].calc(ave_Gs, ZsumB);
	}

	maxGx = profile[0].get_ave_Gx();
	for(int j = 1; j < length_seq; j++){
		if(maxGx < profile[j].get_ave_Gx()){
			maxGx = profile[j].get_ave_Gx(); // find max Gx
		}
	}

	for(int k = 0; k < length_seq; k++){
		if(profile[k].get_ave_Gx() == INFINITE_Gx) {
            profile[k].set_ave_Gx(maxGx); // replace INF Gx with maxGx
        }
	}
}

void SIDD_Base::sum_open()
{
	double sumn=0.0;
	for(int i = 0; i < length_seq; i++)
 	sumn+=profile[i].get_px();
	if (results) {
		if (EnergyType == Z_DNA)
			cout << "Number of Z-DNA bases = " << sumn << endl;
		else if (EnergyType == Cruciform)
			cout << "Number of cruciform bases = " << sumn << endl;
		else
			cout << "Number of melted bases = " << sumn << endl;
	}
}

void SIDD_Base::get_column_header()
{
    if(EnergyType==Z_DNA || EnergyType == Cruciform){
        if(get_showbase()) {
            cout << "Position" << "\tBase" << "\tP(x)"<<endl;
        }
        else {
           cout << "Position" << "\tP(x)"<<endl; 
        }
	}
	else {
        if(get_showbase()) {
             cout << "Position" << "\tBase" << "\tP(x)"<<"\tG(x)"<<endl;
         }
         else {
             cout << "Position" << "\tP(x)"<<"\tG(x)"<<endl;
         }
    }
}

void SIDD_Base::show_profile()
{
	if(MoleculeType == Linear){
		length_seq = length_seq - Len_Gap_Seq;
	}
	get_column_header();
	double Adj_Gx = maxGx - GMAX_ADJ;
	for(int i = 0; i < length_seq; i++){
        if(EnergyType==Z_DNA || EnergyType == Cruciform){
            if(get_showbase()) {
                cout << i + 1 << "\t"<<decode_base(plasmid_seq[i]) <<"\t"<< profile[i].get_px()<<endl;
            }
            else{
                cout << i + 1 <<"\t"<< profile[i].get_px()<<endl;
        	}
        }
        else{
            if(get_showbase()) {
                cout << i + 1 << "\t"<<decode_base(plasmid_seq[i])<<"\t"<< profile[i].get_px()<<"\t"<<profile[i].get_ave_Gx()-Adj_Gx<<endl;
            }
            else {
                cout << i + 1 <<"\t"<< profile[i].get_px()<<"\t"<<profile[i].get_ave_Gx()-Adj_Gx<<"\t"<<endl;
                
            }
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
			cout  << 'C';
			break;
		case 2: 
			cout << 'G';
			break;
		case 3: 
			cout << 'T';
			break;
		default: 
			cout  << 'N';
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


void SIDD_Base::show_parameter()
{
	
	cout << "Sequence Length = " << get_sequence_length() << endl;
	if(EnergyType == Copolymeric)
		cout << "Transition: Melting - copolymeric energetics\n";
	else if(EnergyType == Near_Neighbor)
		cout << "Energetics: Melting - nearest neighbor energetics\n";
	else if (EnergyType == Z_DNA)
		cout << "Transition: Z-DNA\n";
	else if (EnergyType == Cruciform)
        cout << "Transition: Cruciform\n";
	if(MoleculeType == Circular)
		cout << "Molecule type: circular DNA\n";
	else
		cout << "Molecule type: linear DNA\n";
	cout << "Stress Level = " << supdensity << endl;
	cout << "linking difference = " << alpha << " TURNS" << endl;
    cout << "Threshold = " << theta << " kcal/mol" << endl;
	cout << "Temperature = " << Temperature << " K" << endl;
	cout << "Salt Concentration = " << Salt_Conc << " M" << endl;
	cout << "C = " << C << " kcal-bp/rad^2" << endl;
    if (EnergyType == Cruciform) {
        set_Cinitiation();
        cout << "Nucleation energy = " << Ecr << " kcal/mol" << endl;
    }
	else
        cout << "Nucleation energy = " << a << " kcal/mol" << endl;
	if (EnergyType != Z_DNA) {	
	cout << "BAT = " << BAT << " kcal/mol; BGC = " << BGC << " kcal/mol" << endl;
	}
	cout << "MinE = " << min_E << " kcal/mol; Max_E = " << max_E << " kcal/mol" << endl;
}

void SIDD_Base::show_avg_runs()
{
	if(ZsumB != 0.0){
		double runs= (sum_1RB + 2*sum_2RB+3*sum_3RB+4*sum_4RB)/ZsumB;    
		cout << "Average number of runs = " << runs << endl;
	
	}
}

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
