/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Young-Keun Kim
Created         : 01-04-2019
Modified        : 11-27-2023
Language/ver	: C in MSVS2017
Course			: Numerical Programming

Description     : TU System of NonLinear
/------------------------------------------------------------------------------------------*/


#include "stdio.h"
#include "stdlib.h"
#include <math.h>

#include "../../include/myNP_22100252.h"
#include "../../include/myMatrix_22100252.h"





// Goes to myNP.h



// Defined in main code
Matrix myFuncEx1(Matrix X);
Matrix myJacobEx1(Matrix X);
Matrix myFuncEx2(Matrix X);
Matrix myJacobEx2(Matrix X);



void main() {

	/*==========================================================================*/
	/*					Variables declaration & initialization					*/
	/*--------------------------------------------------------------------------*/
	/*   - You can change the variable names									*/
	/*   - However, you must use the specified file name						*/
	/*	   : For each assignment, the file name will be notified on HISNET		*/
	/*==========================================================================*/

	/************      Variables declaration & initialization      ************/
	double loss = 0;
	double n1 = 2;
	double n2 = 3;

	Matrix J = zeros(n1, n1);
	Matrix F = zeros(n1, 1);
	Matrix H = zeros(n1, 1);
	Matrix Z_Q1 = zeros(n1, 1);
	Matrix Z_Q1_initial = zeros(n1, 1);

	Matrix Z_Q2 = zeros(n2, 1);
	Matrix Z_Q2_initial = zeros(n2, 1);

	// Initial condition
	double z0_Q1[2] = { 2.5, 2.0 };
	Z_Q1_initial = arr2Mat(z0_Q1, n1, 1);


	double z0_Q2[3] = { 25.0 / 180.0 * PI, 90.0, 90.0 } ; 

	Z_Q2_initial = arr2Mat(z0_Q2, n2, 1) ;








	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/
	printf("----------------------------------------------------------------------------------------------\n");
	printf("			            System of NonLinear Q1													  \n");
	printf("----------------------------------------------------------------------------------------------\n\r");


	// Áß¿ä
	Z_Q1 = nonlinearSys(myFuncEx1, myJacobEx1, Z_Q1_initial, 0.001);
	printMat(Z_Q1, "Z_Q1");


	printf("----------------------------------------------------------------------------------------------\n");
	printf("			            System of NonLinear Q2													  \n");
	printf("----------------------------------------------------------------------------------------------\n\r");


	Z_Q2 = nonlinearSys(myFuncEx2, myJacobEx2, Z_Q2_initial, 0.00001);
	printMat(Z_Q2, "Z_Q2");





	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/
	freeMat(J);	freeMat(H);	freeMat(F); 	freeMat(Z_Q1);

	system("pause");
}




/*==========================================================================*/
/*						Function Definitions								*/
/*==========================================================================*/

Matrix myFuncEx1(Matrix X)
{

	int n = X.rows;
	Matrix F = zeros(n, 1);
	double x1 = X.at[0][0];
	double x2 = X.at[1][0];

	// [TO-DO] YOUR CODE GOES HERE

	F.at[0][0] = (x2) - 0.5 * (exp(0.5*(x1)) + exp(0.5*(-x1))) ;
	F.at[1][0] = 9 * (x1) * (x1) + 25.0 * x2 * x2 - 225.0 ; 

	return F ;

}


Matrix myJacobEx1(Matrix X)
{
	int n = X.rows;
	Matrix J = zeros(n, n);
	double x1 = X.at[0][0];
	double x2 = X.at[1][0];

	// [TO-DO] YOUR CODE GOES HERE
	J.at[0][0] = - 0.25 * ( exp(0.5 * (x1)) -  exp(0.5*(-x1)) ) ; 
	J.at[0][1] = 1.0 ; 
	J.at[1][0] = 2.0 * 9.0 * (x1) ; 
	J.at[1][1] = 50.0 * (x2) ; 

	return J;
}


Matrix myFuncEx2(Matrix X) {

	int n = X.rows;

	Matrix F = zeros(n, 1); 

	double x1 = X.at[0][0]; 
	double x2 = X.at[1][0]; 
	double x3 = X.at[2][0]; 

	double P0x = 0.0 ;
	double P0y = 100.0 ;
	double P1x = 0.0 ;
	double P1y = -100.0 ;

	double P0x_new = 50.0 ;
	double P0y_new = 186.6025 ;
	double P1x_new = 150.0 ;
	double P1y_new = 13.3975 ;

	// [TO-DO] YOUR CODE GOES HERE
	F.at[0][0] = P0x * cos(x1) - P0y * sin(x1) + (x2) - P0x_new ; 
	F.at[1][0] = P1x * cos(x1) - P1y * sin(x1) + (x2) - P1x_new ;
	F.at[2][0] = P0x * sin(x1) + P0y * cos(x1) + (x3) - P0y_new ; 
	

	return F ;
}

Matrix myJacobEx2(Matrix X) {

	int n = X.rows;

	Matrix J = zeros(n, n) ;
	double x1 = X.at[0][0] ;
	double x2 = X.at[1][0] ;
	double x3 = X.at[2][0] ;

	double P0x = 0.0 ;    
	double P0y = 100.0 ;
	double P1x = 0.0 ;          
	double P1y = -100.0 ;


	J.at[0][0] = -P0x * sin(x1) - P0y * cos(x1) ; 
	J.at[0][1] = 1.0 ; 
	J.at[0][2] = 0.0 ; 

	J.at[1][0] = -P1x * sin(x1) - P1y * cos(x1) ; 
	J.at[1][1] = 1.0 ; 
	J.at[1][2] = 0.0 ; 

	J.at[2][0] = P0x * cos(x1) - P0y * sin(x1) ; 
	J.at[2][1] = 0.0 ; 
	J.at[2][2] = 1.0 ; 

	return J;

}



