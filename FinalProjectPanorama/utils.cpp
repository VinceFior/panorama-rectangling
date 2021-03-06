//
//  utils.cpp
//  FinalProjectPanorama
//
//  Created by Vincent Fiorentini on 4/28/15.
//  Copyright (c) 2015 Vincent Fiorentini. All rights reserved.
//

#include "utils.h"
#include <math.h>
#include <iostream>

// Some methods taken from https://chi3x10.wordpress.com/2008/05/28/calculate-matrix-inversion-in-c/ ,
//  and http://hullooo.blogspot.com/2011/02/matrix-inversion-by-gauss-jordan.html with
//  https://www.youtube.com/watch?v=kMApKEKisKE .

/*
 * The given matrix must be square.
 *
 */
vector<vector<double>> matinv(vector<vector<double>> m, bool useGJ)
{
    int order = (int) m.size();
    
    double **A; // input (m)
    double **Y; // output
    A = new double*[order];
    Y = new double*[order];
    for (int i = 0; i < order; i++) {
        A[i] = new double[order];
        Y[i] = new double[order];
    }
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            A[i][j] = m[i][j];
        }
    }
    if (useGJ) {
        MatrixInversionGJ(A, order, Y);
    } else {
        MatrixInversionMinor(A, order, Y);
    }
    vector<vector<double>> result;
    for (int i = 0; i < order; i++) {
        vector<double> row;
        for (int j = 0; j < order; j++) {
            row.push_back(Y[i][j]);
        }
        result.push_back(row);
    }
    return result;
}

vector<vector<double>> matinvMinor(vector<vector<double>> matrix)
{
    return matinv(matrix, false);
}

vector<vector<double>> matinvGJ(vector<vector<double>> matrix)
{
    return matinv(matrix, true);
}

// Matrix inversion, using minors.
// The result is put in Y
void MatrixInversionMinor(double **A, int order, double **Y)
{
    // get the determinant of a
    double det = 1.0/CalcDeterminant(A,order);
    
    // memory allocation
    double *temp = new double[(order-1)*(order-1)];
    double **minor = new double*[order-1];
    for(int i=0;i<order-1;i++)
        minor[i] = temp+(i*(order-1));
    
    for(int j=0;j<order;j++)
    {
        for(int i=0;i<order;i++)
        {
            // get the co-factor (matrix) of A(j,i)
            GetMinor(A,minor,j,i,order);
            Y[i][j] = det*CalcDeterminant(minor,order-1);
            if( (i+j)%2 == 1)
                Y[i][j] = -Y[i][j];
        }
    }
    
    // release memory
    //delete [] minor[0];
    delete [] temp;
    delete [] minor;
}

// calculate the cofactor of element (row,col)
int GetMinor(double **src, double **dest, int row, int col, int order)
{
    // indicate which col and row is being copied to dest
    int colCount=0,rowCount=0;
    
    for(int i = 0; i < order; i++ )
    {
        if( i != row )
        {
            colCount = 0;
            for(int j = 0; j < order; j++ )
            {
                // when j is not the element
                if( j != col )
                {
                    dest[rowCount][colCount] = src[i][j];
                    colCount++;
                }
            }
            rowCount++;
        }
    }
    
    return 1;
}

// Calculate the determinant recursively.
double CalcDeterminant( double **mat, int order)
{
    // order must be >= 0
    // stop the recursion when matrix is a single element
    if( order == 1 )
        return mat[0][0];
    
    // the determinant value
    double det = 0;
    
    // allocate the cofactor matrix
    double **minor;
    minor = new double*[order-1];
    for(int i=0;i<order-1;i++)
        minor[i] = new double[order-1];
    
    for(int i = 0; i < order; i++ )
    {
        // get minor of element (0,i)
        GetMinor( mat, minor, 0, i , order);
        // the recusion is here!
        
        det += (i%2==1?-1.0:1.0) * mat[0][i] * CalcDeterminant(minor,order-1);
        //det += pow( -1.0, i ) * mat[0][i] * CalcDeterminant( minor,order-1 );
    }
    
    // release memory
    for(int i=0;i<order-1;i++)
        delete [] minor[i];
    delete [] minor;
    
    return det;
}

// Matrix inversion, using Gauss-Jordan elimination.
// The result is put in Y.
void MatrixInversionGJ(double **A, int order, double **Y)
{
    double **augmentedmatrix;
    augmentedmatrix = new double*[order];
    for (int i = 0; i < order; i++) {
        augmentedmatrix[i] = new double[2 * order];
    }
    
    // store the augmented matrix as a matrix of dimensions (order)x(2*order) in a 2D array
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            augmentedmatrix[i][j] = A[i][j];
        }
    }
    
    // augment with the identity matrix
    for (int i = 0; i < order; i++) {
        for (int j = order; j < 2 * order; j++) {
            if (i == (j % order)) {
                augmentedmatrix[i][j] = 1;
            } else {
                augmentedmatrix[i][j] = 0;
            }
        }
    }
    
//    cerr << "The initial augmented matrix is: " << endl;
//    for(int i = 0; i < order; i++) {
//        for(int j = 0; j < 2 * order; j++) {
//            cerr << " " << augmentedmatrix[i][j];
//        }
//        cerr << endl;
//    }
    
    // convert to upper triangular
    for (int i = 0; i < order - 1; i++) {
        // i is the row (so like y, from the top)
        for (int k = 1; k < order; k++) {
            // k is the row to switch with
            if (augmentedmatrix[i][i] == 0) {
                // swap rows
                for (int l = 0; l < 2 * order; l++) {
                    double tmp = augmentedmatrix[i][l];
                    augmentedmatrix[i][l] = augmentedmatrix[k][l];
                    augmentedmatrix[k][l] = tmp;
                }
            }
        }
        for (int j = i + 1; j < order; j++) {
            // subtract other row
            double factor = augmentedmatrix[j][i] / augmentedmatrix[i][i];
            for (int l = 0; l < 2 * order; l++) {
                augmentedmatrix[j][l] -= augmentedmatrix[i][l] * factor;
            }
        }
    }

    // convert to diagonal
    for (int j = 0; j < order; j++) {
        // j is the column
        for (int i = 0; i < j; i++) {
            // i is the row
            if (augmentedmatrix[i][j] != 0) {
                // subtract from this row a multiple of a row beneath it
                double factor = -1 * augmentedmatrix[i][j] / augmentedmatrix[j][j];
                for (int l = 0; l < 2 * order; l++) {
                    augmentedmatrix[i][l] += factor * augmentedmatrix[j][l];
                }
            }
        }
    }
    
    // normalize
    for (int i = 0; i < order; i++) {
        // i is the row
        double factor = 1 / augmentedmatrix[i][i];
        for (int l = 0; l < 2 * order; l++) {
            augmentedmatrix[i][l] *= factor;
        }
    }
    
//    cerr << "After Gauss-Jordan elimination, the augmented matrix is: " << endl;
//    for(int i = 0; i < order; i++) {
//        for(int j = 0; j < 2 * order; j++) {
//            cerr << " " << augmentedmatrix[i][j];
//        }
//        cerr << endl;
//    }

//    cerr << "The inverse of the given matrix is: " << endl;
//    for(int i = 0; i < order; i++) {
//        for(int j = order; j < 2 * order; j++) {
//            cerr << " " << augmentedmatrix[i][j];
//        }
//        cerr << endl;
//    }
    
    // move the right side of the augmented matrix into the output array Y
    for (int i = 0; i < order; i++) {
        for (int j = order; j < 2 * order; j++) {
            Y[i][j - order] = augmentedmatrix[i][j];
        }
    }
    
//    cerr << "The output is: " << endl;
//    for (int i = 0; i < order; i++) {
//        for (int j = 0; j < order; j++) {
//            cerr << " " << Y[i][j];
//        }
//        cerr << endl;
//    }
}