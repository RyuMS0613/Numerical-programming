/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Min-Seo Ryu - Handong Global University

Author          : Min-Seo Ryu
Created         : 09-07-2024
Modified        : 09-11-2024
Language/ver    : C /  MSVS2017
Course          : Numerical method 2024-2
Description     : Assignment of Newton Rhapson  (Solving Non-linear Eqation)
/------------------------------------------------------------------------------------------*/



#include "../../include/myNP_22100252.h"

/* function f(x) of the problem */
double func(double x);

/* first derivative of function f(x) of the problem */
double dfunc(double x);


void main() {

	/*==========================================================================*/
	/*               Assignment -     Newton Rhapson                            */
	/*==========================================================================*/

	/************      Variables declaration & initialization      ************/
	double tol = 1e-6;
	double sol_nr = 0;

	double a = 27.00 ;
	double b = 29.00 ;

	printf("------------------------------------------------------------------------------------\n");
	printf("         Newton-Raphson Method Results             \n");
	printf("------------------------------------------------------------------------------------\n");

	printf("Secant Method Result:\n");
	printf("Final Solution: %f \t", bisection(func, a, b, tol));
	printf("\n");


	double x0 = 28.00;
	double x0_1 = 27.00;



	printf("------------------------------------------------------------------------------------\n");
	printf("         Newton-Raphson Method Results             \n");
	printf("------------------------------------------------------------------------------------\n");


	sol_nr = newtonRaphson(func, dfunc, x0, tol); 

	printf("Newton-Raphson Method Result:\n") ;
	printf("Final Solution: %f \t", sol_nr);
	printf("\n");


	printf("------------------------------------------------------------------------------------\n");
	printf("         Secant Method Results             \n");
	printf("------------------------------------------------------------------------------------\n");

	sol_nr = secant(func, x0, x0_1, tol);

	printf("Secant Method Result:\n"); 
	printf("Final Solution: %f \t", sol_nr); 
	printf("\n"); 


	system("pause");
}




/*==========================================================================*/
/*                    Function Definitions		   		*/
/*==========================================================================*/


double func(double x)
{
	double F = 0;

	F = 20 * tan(x) - (9.8 * power(20, 2)) / (2 * power(17, 2) * power(cos(x), 2)) - 2;


	return F;
}

/* first derivative of function f(x) of the problem */
double dfunc(double x)
{

	double dF = 0;

	dF = 20 / power(cos(x), 2) - (9.8 * power(20, 2) * sin(x)) / (power(17, 2) * power(cos(x), 3));
	
	return dF;
}
