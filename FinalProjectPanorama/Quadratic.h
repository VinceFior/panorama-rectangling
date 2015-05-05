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
#include <vector>

using namespace std;

struct QuadraticTerm
{
    // -1 for either index indicates the term is only in one variable (degree 1)
    int varIndex1;
    int varIndex2;
    double coefficientValue;
};

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
    void incrementCoeffForVar(int varIndex, double value);
    void incrementCoeffForVars(int varIndex1, int varIndex2, double value);
    double **derivatives();
    void addQuadratic(Quadratic quad, double weight);
    void scaleCoefficients(double scaleFactor);
    void roundCoefficients(double epsilon);
    void printEquation();
private:
    vector<QuadraticTerm> coefficients;
    int numVariables; // the number of variables in the input domain
    size_t quadraticIndexForInputIndex(int inputIndex);
    size_t quadraticIndexForInputIndices(int inputIndex1, int inputIndex2);
    double *derivativeWithRespectToInputIndex(int inputIndex);
    string varIndexToStr(int varIndex);
};

#endif /* defined(__FinalProjectPanorama__Quadratic__) */
