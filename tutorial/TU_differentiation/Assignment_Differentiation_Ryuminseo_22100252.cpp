/*------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author          : Min-Seo Ryu
Created         : 09-16-2024
Modified        : 09-17-2024
Language/ver	: C in MSVS2017
Course		: Numerical Programming

Description      : [Assignment] Numerical Differntiation.cpp
-------------------------------------------------------------------------------*/


#include "../../include/myNP_22100252.h"

double myFunc(const double x);


int main(int argc, char* argv[])
{

	/*==========================================================================*/
	/*   Part 1 -     Differentiation from discrete dataset points              */
	/*==========================================================================*/

	printf("\n**************************************************");
	printf("\n|                     PART 1.                    |");
	printf("\n**************************************************\n");

	/************      Variables declaration & initialization      ************/
	int m1 = 12;
	double t1[12] = { -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5 };

	double x[12] = { -3.632, - 0.3935, 1.0, 0.6487, - 1.282, - 4.518, - 8.611, - 12.82, - 15.91, - 15.88, - 9.402, 9.017 };
	double  dxdt1[12] = { 0 };

	/************      Solve  &	Show Output	   ************/
	// Differentiation from discrete dataset points

	// [YOUR CODE GOES HERE]
	gradient1D(t1, x, dxdt1, m1);
	printVec(dxdt1, m1);



	/*==========================================================================*/
	/*   Part 2 -     Differentiation from a function                           */
	/*==========================================================================*/


	printf("\n**************************************************");
	printf("\n|                     PART 2.                    |");
	printf("\n**************************************************\n");

	/************      Variables declaration & initialization      ************/
	double xin = 2;
	int m2 = 21;
	double dydx2[21] = { 0 };  // m=12 points
	double t2[21] = { 0 };

	for (int i = 0; i < 21; i++) {
		t2[i] = i*0.2;
	}

	// User defined function F(x)
	double y = myFunc(xin);
	printf("\n y=myFun(xin) = %f \n\n", y);

	/************      Solve  &	Show Output	   ************/
	// Estimate differentiation from the user defined function 

	// [YOUR CODE GOES HERE]
	gradientFunc(myFunc, t2, dydx2, m2);
	printVec(dydx2, m2);


	/*==========================================================================*/
	/*   Part 3 -     Second Differentiation from a function                           */
	/*==========================================================================*/

	printf("\n**************************************************");
	printf("\n|                     PART 3.                    |");
	printf("\n**************************************************\n");

	int m3 =21;
	double dydt3[21] = { 0 };

	double t3[21] = { 0 };
	for (int i = 0; i < 21; i++) {
		t3[i] = i * 0.2;
	}

	double x3[21] = { 0 };
	for (int i = 0; i < 21; i++) {
		x3[i] = myFunc(t3[i]);
	}


	acceleration(t3, x3, dydt3, m3);
	printVec(dydt3, m3);

	system("pause");
	return 0;
}


// User defined function:  example  y=x*x
// Modify to User Function
double myFunc(const double x) {
	return  power(x,3);
}