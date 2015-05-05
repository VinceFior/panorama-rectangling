//
//  Quadratic.cpp
//  FinalProjectPanorama
//
//  Created by Vincent Fiorentini on 5/4/15.
//  Copyright (c) 2015 Vincent Fiorentini. All rights reserved.
//

#include "Quadratic.h"
#include <iostream>
#include <math.h>

Quadratic::Quadratic(int numVariables)
{
    this->numVariables = numVariables;
}

Quadratic::~Quadratic()
{
    // nothing to do
}

double Quadratic::getCoeffForVar(int varIndex)
{
    for (int i = 0; i < coefficients.size(); i++) {
        QuadraticTerm term = coefficients[i];
        if ((term.varIndex1 == varIndex && term.varIndex2 == -1) ||
            (term.varIndex2 == varIndex && term.varIndex1 == -1)) {
            return term.coefficientValue;
        }
    }
    return 0;
}

double Quadratic::getCoeffForVars(int varIndex1, int varIndex2)
{
    for (int i = 0; i < coefficients.size(); i++) {
        QuadraticTerm term = coefficients[i];
        if ((term.varIndex1 == varIndex1 && term.varIndex2 == varIndex2) ||
            (term.varIndex1 == varIndex2 && term.varIndex2 == varIndex1)) {
            return term.coefficientValue;
        }
    }
    return 0;
}

void Quadratic::setCoeffForVar(int varIndex, double value)
{
    for (int i = 0; i < coefficients.size(); i++) {
        QuadraticTerm& term = coefficients[i];
        if ((term.varIndex1 == varIndex && term.varIndex2 == -1) ||
            (term.varIndex2 == varIndex && term.varIndex1 == -1)) {
            term.coefficientValue = value;
            return;
        }
    }
    QuadraticTerm term;
    term.varIndex1 = varIndex;
    term.varIndex2 = -1;
    term.coefficientValue = value;
    coefficients.push_back(term);
}

void Quadratic::setCoeffForVars(int varIndex1, int varIndex2, double value)
{
    for (int i = 0; i < coefficients.size(); i++) {
        QuadraticTerm& term = coefficients[i];
        if ((term.varIndex1 == varIndex1 && term.varIndex2 == varIndex2) ||
            (term.varIndex1 == varIndex2 && term.varIndex1 == varIndex1)) {
            term.coefficientValue = value;
            return;
        }
    }
    QuadraticTerm term;
    term.varIndex1 = varIndex1;
    term.varIndex2 = varIndex2;
    term.coefficientValue = value;
    coefficients.push_back(term);
}

void Quadratic::incrementCoeffForVar(int varIndex, double value)
{
    setCoeffForVar(varIndex, value + getCoeffForVar(varIndex));
}

void Quadratic::incrementCoeffForVars(int varIndex1, int varIndex2, double value)
{
    setCoeffForVars(varIndex1, varIndex2, value + getCoeffForVars(varIndex1, varIndex2));
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
 * Adds the given quadratic to the this quadratic.
 */
void Quadratic::addQuadratic(Quadratic quad, double weight)
{
    vector<QuadraticTerm> newCoefficients;
    // add the terms in both equations
    for (int i = 0; i < coefficients.size(); i++) {
        QuadraticTerm term = coefficients[i];
        bool quadHasTerm = false;
        for (int j = 0; j < quad.coefficients.size(); j++) {
            QuadraticTerm quadTerm = quad.coefficients[j];
            if ((quadTerm.varIndex1 == term.varIndex1 && quadTerm.varIndex2 == term.varIndex2) ||
                (quadTerm.varIndex1 == term.varIndex2 && quadTerm.varIndex2 == term.varIndex1)) {
                QuadraticTerm newTerm;
                newTerm.varIndex1 = term.varIndex1;
                newTerm.varIndex2 = term.varIndex2;
                newTerm.coefficientValue = term.coefficientValue + weight * quadTerm.coefficientValue;
                newCoefficients.push_back(newTerm);
                quadHasTerm = true;
            }
        }
        // add the terms in only this equation
        if (!quadHasTerm) {
            QuadraticTerm newTerm;
            newTerm.varIndex1 = term.varIndex1;
            newTerm.varIndex2 = term.varIndex2;
            newTerm.coefficientValue = term.coefficientValue;
            newCoefficients.push_back(newTerm);
        }
    }
    // add the terms in only the other quad equation
    for (int i = 0; i < quad.coefficients.size(); i++) {
        QuadraticTerm quadTerm = quad.coefficients[i];
        bool hasQuadTerm = false;
        for (int j = 0; j < coefficients.size(); j++) {
            QuadraticTerm term = coefficients[j];
            if ((quadTerm.varIndex1 == term.varIndex1 && quadTerm.varIndex2 == term.varIndex2) ||
                (quadTerm.varIndex1 == term.varIndex2 && quadTerm.varIndex2 == term.varIndex1)) {
                hasQuadTerm = true;
            }
        }
        // not in this equation, only in the other quad equation
        if (!hasQuadTerm) {
            QuadraticTerm newTerm;
            newTerm.varIndex1 = quadTerm.varIndex1;
            newTerm.varIndex2 = quadTerm.varIndex2;
            newTerm.coefficientValue = weight * quadTerm.coefficientValue;
            newCoefficients.push_back(newTerm);
        }
    }
    
    // set our coefficients to the new coefficients
    coefficients = newCoefficients;
}

/*
 * Multiplies all coefficients by the given scale factor.
 */
void Quadratic::scaleCoefficients(double scaleFactor)
{
    for (int i = 0; i < coefficients.size(); i++) {
        QuadraticTerm& term = coefficients[i];
        term.coefficientValue *= scaleFactor;
    }
}

/*
 * Rounds all coefficients of absolute value less than epsilon to 0.
 */
void Quadratic::roundCoefficients(double epsilon)
{
    for (int i = 0; i < coefficients.size(); i++) {
        QuadraticTerm &term = coefficients[i];
        if (fabs(term.coefficientValue) < epsilon) {
            term.coefficientValue = 0;
        }
    }
}

/*
 * Unfinished convenience method.
 */
string Quadratic::varIndexToStr(int varIndex)
{
    if (varIndex == 0) {
        return "x0";
    } else if (varIndex == 1) {
        return "y0";
    } else if (varIndex == 2) {
        return "x1";
    } else if (varIndex == 3) {
        return "y1";
    } else if (varIndex == 4) {
        return "x2";
    } else if (varIndex == 5) {
        return "y2";
    } else if (varIndex == 6) {
        return "x3";
    } else if (varIndex == 7) {
        return "y3";
    } else {
        return to_string(varIndex);
    }
}

/*
 * Prints out the coefficients of the variables in the equation (to cerr).
 * Ex., the equation 7a^2 + 2ab - 3b^2 + 5a is 7:0*0 2:0*1 -3:1*1 5:0.
 */
void Quadratic::printEquation()
{
    for (int i = 0; i < coefficients.size(); i++) {
        QuadraticTerm term = coefficients[i];
        if (term.varIndex1 != -1 && term.varIndex2 != -1) {
            cerr << term.coefficientValue << ":" << varIndexToStr(term.varIndex1) << "*" << varIndexToStr(term.varIndex2) << " ";
        } else if (term.varIndex1 != -1 && term.varIndex2 == -1) {
            cerr << term.coefficientValue << ":" << varIndexToStr(term.varIndex1) << " ";
        } else if (term.varIndex2 != -1 && term.varIndex1 == -1) {
            cerr << term.coefficientValue << ":" << varIndexToStr(term.varIndex2) << " ";
        } else {
            // this case should never happen (because we don't include constant terms)
            cerr << term.coefficientValue << ":: ";
        }
    }
    cerr << endl;
}
