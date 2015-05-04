//
//  Quadratic.h
//  FinalProjectPanorama
//
//  Created by Vincent Fiorentini on 5/4/15.
//  Copyright (c) 2015 Vincent Fiorentini. All rights reserved.
//

#ifndef __FinalProjectPanorama__Quadratic__
#define __FinalProjectPanorama__Quadratic__

#include <stdio.h>

// Note: This representation of a quadratic equation does not allow for constant terms.
// Also, this class is hard-coded to represent quadratic equations with double coefficients.
class Quadratic
{
public:
    Quadratic (int numVariables);
    ~Quadratic ();
    double getCoeffForVar(int varIndex);
    double getCoeffForVars(int varIndex1, int varIndex2);
    void setCoeffForVar(int varIndex, double value);
    void setCoeffForVars(int varIndex1, int varIndex2, double value);
    double **derivatives();
    void addQuadratic(Quadratic quad, double weight);
    void printEquation();
private:
    double *equation; // the coefficients for each term
    int numVariables; // the number of variables in the input domain
    size_t numTerms; // the number of terms in the quadratic; no constants
    size_t quadraticIndexForInputIndex(int inputIndex);
    size_t quadraticIndexForInputIndices(int inputIndex1, int inputIndex2);
    double *derivativeWithRespectToInputIndex(int inputIndex);
};

#endif /* defined(__FinalProjectPanorama__Quadratic__) */
