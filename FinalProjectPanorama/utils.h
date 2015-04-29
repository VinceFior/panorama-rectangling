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
void MatrixInversion(float **A, int order, float **Y);
int GetMinor(float **src, float **dest, int row, int col, int order);
double CalcDeterminant( float **mat, int order);

#endif /* defined(__FinalProjectPanorama__utils__) */
