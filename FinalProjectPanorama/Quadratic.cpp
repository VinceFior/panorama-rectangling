//
//  Quadratic.cpp
//  FinalProjectPanorama
//
//  Created by Vincent Fiorentini on 5/4/15.
//  Copyright (c) 2015 Vincent Fiorentini. All rights reserved.
//

#include "Quadratic.h"
#include <iostream>

/*
 * Computes n!
 */
size_t factorial(size_t n)
{
    if (n == 0) {
        return 1;
    } else {
        return n * factorial(n - 1);
    }
}

Quadratic::Quadratic(int numVariables)
{
    // the terms are (not in order): each variable squared, each two choices of variables, each one variable
    // numTerms = numVariables + (numVariables choose 2) + numVariables
    //  = n!/(2(n-2)!) + 2 * numVariables
    this->numTerms = factorial(numVariables) / (2 * factorial(numVariables - 2)) + 2 * numVariables;
    this->numVariables = numVariables;
    this->equation = new double[numTerms];
    for (int i = 0; i < numTerms; i++) {
        equation[i] = 0;
    }
}

Quadratic::~Quadratic()
{
    // nothing to do
}

/* The order of the terms is: var_i * var_j for i=0 to i=numVariables (for j=0 to j=numVariables),
 * then var_i for i=0 to i=numVariables.
 * Ex., for variables a,b,c, the order is: aa, ab, ac, bb, bc, cc, a, b, c.
 *   In this example, input index 2 (b) yields quadratic index 7 (b).
 */
size_t Quadratic::quadraticIndexForInputIndex(int inputIndex)
{
    // the first term that is a single variable is after numVariables + (numVariables choose 2)
    size_t firstIndex = factorial(numVariables) / (2 * factorial(numVariables - 2)) + numVariables;
    return firstIndex + inputIndex;
}

/* We use the upper-right triangle of the matrix [a, b, c] * [a, b, c]^T
 *
 *   aa ab ac
 *   ba bb bc
 *   ca cb cc
 * Ex., the index for bc is 4.
 */
size_t Quadratic::quadraticIndexForInputIndices(int inputIndex1, int inputIndex2)
{
    int smallTermIndex;
    int largeTermIndex;
    if (inputIndex1 > inputIndex2) {
        largeTermIndex = inputIndex1;
        smallTermIndex = inputIndex2;
    } else {
        largeTermIndex = inputIndex2;
        smallTermIndex = inputIndex1;
    }
    size_t index = 0;
    for (int i = 0; i < smallTermIndex; i++) {
        index += numVariables - i;
    }
    index += (largeTermIndex - smallTermIndex);
    return index;
}

// varIndex is the index of the variable in the input list.
// Ex., for input variables a,b,c, the index of variable b is 1.
double Quadratic::getCoeffForVar(int varIndex)
{
    return equation[quadraticIndexForInputIndex(varIndex)];
}

double Quadratic::getCoeffForVars(int varIndex1, int varIndex2)
{
    return equation[quadraticIndexForInputIndices(varIndex1, varIndex2)];
}

void Quadratic::setCoeffForVar(int varIndex, double value)
{
    equation[quadraticIndexForInputIndex(varIndex)] = value;
}

void Quadratic::setCoeffForVars(int varIndex1, int varIndex2, double value)
{
    equation[quadraticIndexForInputIndices(varIndex1, varIndex2)] = value;
}

/*
 * Each of the numVariables rows is a linear equation of numVariables+1 elements.
 * (There is one row (equation) for each variable, and each equation has numVariables + 1 terms,
 *  one for each variable and one for a constant.)
 */
double **Quadratic::derivatives()
{
    double **derivatives = new double* [numVariables];
    for (int i = 0; i < numVariables; i++) {
        derivatives[i] = derivativeWithRespectToInputIndex(i);
    }
    return derivatives;
}

/*
 * Returns a linear equation (specifically, the coefficients of the input variables, in order).
 * Ex., the output order is a,b,c,const.
 */
double *Quadratic::derivativeWithRespectToInputIndex(int inputIndex)
{
    double *derivative = new double[numVariables + 1];
    // all terms are set to 0 by default
    for (int i = 0; i < numVariables + 1; i++) {
        derivative[i] = 0;
    }
    // squared term of interest
    derivative[inputIndex] = 2 * getCoeffForVars(inputIndex, inputIndex);
    // non-squared quadratic terms
    for (int i = 0; i < numVariables; i++) {
        for (int j = i + 1; j < numVariables; j++) {
            if (i == inputIndex) {
                derivative[j] = getCoeffForVars(i, j);
            } else if (j == inputIndex) {
                derivative[i] = getCoeffForVars(i, j);
            }
        }
    }
    // constant term
    derivative[numVariables] = getCoeffForVar(inputIndex);
    
    return derivative;
}

/*
 * Adds the given quadratic to the this quadratic. Assumes the quadratics are of the same number of
 * variables.
 */
void Quadratic::addQuadratic(Quadratic quad, double weight)
{
    for (int i = 0; i < numTerms; i++) {
        equation[i] += weight * quad.equation[i];
    }
}

/*
 * Prints out the coefficients of the variables in the equation (to cerr).
 */
void Quadratic::printEquation()
{
    for (int i = 0; i < numTerms; i++) {
        std::cerr << equation[i] << " ";
    }
    std::cerr << std::endl;
}
