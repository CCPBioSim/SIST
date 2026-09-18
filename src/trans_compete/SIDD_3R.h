// SIDD_3R.h: interface for the SIDD_3R class.
//
// Author: Chengpeng Bi
// Last Modified: August 1, 2002
//
// This is a C++ implementation of the SIDD algorithm designed by Dr. Benham
//
// This class is designed for three runs derived from the class of two runs
//
// Detailed documentation to come
//
// UC Davis Genome Center
///////////////////////////////////////////////////////////////////////////////////////

#ifndef SIDD_3R_H_
#define SIDD_3R_H_

#include "SIDD_2R.h"

// Definition of the SIDD_3R class, inheriting from SIDD_2R
class SIDD_3R : public SIDD_2R
{
protected:
    // Variables specific to SIDD_3R class
    bool flag_minE3;
    long StatesInThreeRuns;

public:
    // Constructor for SIDD_3R class
    SIDD_3R();

    // Destructor for SIDD_3R class
    virtual ~SIDD_3R();

    // Function to search for low-energy three-run states
    bool Search_Low3RE();
};

#endif // !defined(AFX_SIDD_3R_H__B5CD0095_10DE_4C22_A817_256ABBAD9425__INCLUDED_)
