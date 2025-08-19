/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 27-09-2024
Modified         : 27-09-2024
Language/ver     : C++ in MSVS2019

Description      : Assignment Numerical Integral
-------------------------------------------------------------------------------*/

#include "stdio.h"
#include "stdlib.h"

#include "../../include/myNP_22100252.h"

// Sample Code: Integration rectangular method 
double IntegrateRect(double x[], double y[], int m);

// You need to create myFunc() in this main source file
double myFunc(const double x);


// You need to create the followins in  in myNP.h. myNP.cpp
// double trapz (double x[ ], double y[ ], int m);
// double simpson13(double x[ ], double y[ ], int m);
// double integral(double func(const double x), double a, double b, int n);  


int main(int argc, char* argv[])
{

	printf("\n**************************************************");
	printf("\n        PART 1. Integration from Datasets         ");
	printf("\n**************************************************\n");

	/************      Variables declaration & initialization      ************/
	double x1[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 };
	double y1[] = { 0, 3, 8, 20, 33, 42, 40, 48, 60, 12, 8, 4, 3 };
	int M1 = sizeof(x1) / sizeof(x1[0]);


	double x2[] = { -3.0 , - 2.25 , -1.50 , -0.75  , 0.0 , 0.75   , 1.5  , 2.25   , 3.0 };
	double y2[] = { 0.0  , 2.1875 ,  3.75 , 4.6875 , 5.0 , 4.6875 , 3.75 , 2.1875 , 0.0 };
	int M2 = sizeof(x2) / sizeof(x2[0]);

	/************      Solve  &	Show Output	   ************/

	// Exercise 1. Trapezoid

	double I_trapz = 0;
	I_trapz = trapz(x1, y1, M1);

	printf("I_trapz = %f\n\n", I_trapz);


	// Exercise 2. Simpson

	double I_simpson13 = 0;
	I_simpson13 = simpson13(x2, y2, M2);

	printf("I_simpson13  = %f\n\n", I_simpson13);




	printf("\n**************************************************");
	printf("\n        PART 2. Integration from a Function       ");
	printf("\n**************************************************\n");

	// Exercise 3. Integral

	double I_function = 0;
	int N = 14;

	I_function = integral(myFunc,-1,1,N);
	printf("I_function  = %f\n\n", I_function);

	system("pause");
	return 0;
}



double myFunc(const double x) {

	return sqrt(1-x*x);
}