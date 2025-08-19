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







int main(int argc, char* argv[])
{

	/*==========================================================================*/
	/*   Part 1 -     Line Spine         */
	/*==========================================================================*/

	printf("\n**************************************************");
	printf("\n|                     PART 1.                    |");
	printf("\n**************************************************\n");


	int m = 21;
	double input = 5.2;

	double x[21] = { 0 };
	for (int i = 0; i < 21; i++) {
		x[i] = i * 1;
	}

	double y[21] = { 0 };
	for (int i = 0; i < 21; i++) {
		y[i] = i * 10;
	}

	//printf("result of lagrange is %f \n", Linespline_lagrange(x, y, input, m));

	//printf("result of Newtonian is %f \n", Linespline_Newtonian(x, y, input, m));

	//printf("result of Newtonian is %f \n", Linespline_Standard(x, y, input, m));

	//printf("result of lagrange is %f \n", lagrange(x,y,m,3,input));

	system("pause");

	double f[10][10] = {0.0};
	double y[10][10] = {0.0};


	for (int i = 0; i < 10; i++) {

		for (int j = 0; j < 10; j++) {
			f[i][j] = i+i*j;
		}
	}
	

	for (int i = 0; i < 10; i++) {

		for (int j = 0; j < 10; j++) {
			printf("%f ", f[i][j]);
		}
		printf("\n");
	}



	return 0;
}

