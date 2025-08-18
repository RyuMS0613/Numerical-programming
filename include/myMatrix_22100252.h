/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 14-10-2024
Language/ver     : C++ in MSVS2019

Description      : myMatrix.h
----------------------------------------------------------------*/

#ifndef		_MY_MATRIX_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_MATRIX_H

#include <iostream>
#include <string>
#include <fstream>

typedef struct { 
	double** at;
	int rows, cols;
}Matrix;


// Create Matrix with specified size
extern	Matrix	createMat(int _rows, int _cols);

// Free a memory allocated matrix
extern	void	freeMat(Matrix _A);

// Create a matrix from a text file
extern	Matrix	txt2Mat(std::string _filePath, std::string _fileName);

//// Print matrix
extern	void	printMat(Matrix _A, const char* _name);

// Matrix addition
extern	Matrix	addMat(Matrix _A, Matrix _B);


// initialization of Matrix elements
extern	void	initMat(Matrix _A, double _val);



//////////////////////////////////////////////////////////////////
/*							Tutorial							*/
//////////////////////////////////////////////////////////////////


// Create matrix of all zeros
extern Matrix zeros(int _rows, int _cols);

// Create matrix of all ones
extern	Matrix	ones(int _rows, int _cols);

// Create identity matrix
extern	Matrix	eye(int _rows, int _cols);


// Matrix subtraction
extern	Matrix	subMat(Matrix _A, Matrix _B);

// Multiply  matrix A and matrix B
extern	Matrix	multMat(Matrix _A, Matrix _B);

extern	void multMat_void(Matrix _A, Matrix _B, Matrix _AB);

// Multiply  matrix A with a scalar k
extern	Matrix	smultMat(Matrix _A, double _k);

// Create Transpose matrix
extern	Matrix	transpose(Matrix _A);

// Copy matrix
void copyMat(Matrix _A, Matrix _B);




////    MATRIX    /////

// Guass Elimination
void gaussElim_Basic(Matrix A, Matrix b, Matrix U, Matrix d);

// Guass Elimination with pivoting
void gaussElim(Matrix A, Matrix b, Matrix U, Matrix d, Matrix P);

// backsub 
void backsub(Matrix U, Matrix d, Matrix x);



// LU decomposition 
void LUdecomp_pivoting(Matrix A, Matrix U, Matrix L,Matrix P) ;

// LU decomposition 
void LUdecomp(Matrix A, Matrix U, Matrix L);

// LU solve 
void LUsolve(Matrix L, Matrix U, Matrix P ,Matrix b, Matrix x) ; 

// forwardsub 
void fwdsub(Matrix L, Matrix d, Matrix x) ; 

// solve Ax = b (with pivoting) 
Matrix solve_system(Matrix A, Matrix b);

Matrix solve_system_pivoting(Matrix A, Matrix b); 

// invert Matrix 
void invMat(Matrix A, Matrix Ainv) ; 

// LUdecomp with scaled partial pivoting  
void scaled_LUdecomp(Matrix A, Matrix U, Matrix L, Matrix P);



// Eigenvalue && Eigenvector


Matrix eigval(Matrix A); 

Matrix eigvec(Matrix A);

void eig(Matrix A, Matrix V, Matrix D);




void QRdecomp(Matrix A, Matrix Q, Matrix R);

double norm_2(Matrix A);



Matrix	linearFit_mat(Matrix _vecX, Matrix _vecY); 

Matrix	arr2Mat(double* _1Darray, int _rows, int _cols); 

Matrix	polyFit_mat(Matrix _vecX, Matrix _vecY, int n);

Matrix	polyFit_mat_count(Matrix _vecX, Matrix _vecY, int n);


Matrix expFit(Matrix _vecX, Matrix _vecY);


Matrix nonlinearSys(Matrix Funcs(Matrix _Z), Matrix Jacob(Matrix _Z), Matrix _Z0, double tol);



#endif