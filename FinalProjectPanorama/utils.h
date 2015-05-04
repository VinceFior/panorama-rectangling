//
//  utils.h
//  FinalProjectPanorama
//
//  Created by Vincent Fiorentini on 4/28/15.
//  Copyright (c) 2015 Vincent Fiorentini. All rights reserved.
//

#ifndef __FinalProjectPanorama__utils__
#define __FinalProjectPanorama__utils__

#include <stdio.h>
#include <vector>
using namespace std;

vector<vector<double>> matinv(vector<vector<double>> m);
void MatrixInversion(double **A, int order, double **Y);
int GetMinor(double **src, double **dest, int row, int col, int order);
double CalcDeterminant( double **mat, int order);

#endif /* defined(__FinalProjectPanorama__utils__) */
