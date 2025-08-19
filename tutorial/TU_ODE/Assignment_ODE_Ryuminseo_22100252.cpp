/*------------------------------------------------------------------------------------------\
@ Numerical Programming  by Young-Keun Kim - Handong Global University

Author          : Young-Keun Kim
Created         : 01-04-2019
Modified        : 01-04-2024
Language/ver	: C in MSVS2017
Course			: Numerical Programming

Description     : Assignment 8 ODE IVP
/------------------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------------------\
@ Numerical Programming  by Young-Keun Kim - Handong Global University

Author          : Young-Keun Kim
Created         : 01-04-2019
Modified        : 01-04-2024
Language/ver	: C in MSVS2017
Course			: Numerical Programming

Description     : Assignment 8 ODE IVP - SOLUTION
/------------------------------------------------------------------------------------------*/

#define ASGN	8		// enter your assignment number
#define EVAL	0		// set 0


/*------------------------------------------------------------------------------------------*/
/*			MODIFY Path and name of Headerfiles										 		*/
/*------------------------------------------------------------------------------------------*/
#include "../../include/myNP_22100252.h"
#include "../../include/myMatrix_22100252.h"



/*------------------------------------------------------------------------------------------*/
/*			MOVE the following  #define to  myNP.h									 		*/
/*------------------------------------------------------------------------------------------*/
#define PI 3.14159265368979323846264338327950288412



/*------------------------------------------------------------------------------------------*/
/*			DEFINE ode model functions here											 		*/
/*------------------------------------------------------------------------------------------*/
// Problem1: RC circuit
double odeFunc_rc(const double t, const double v);

// Problem2: m-c-k system
void odeFunc_mck(double dYdt[], const double t, const double Y[]);



/*------------------------------------------------------------------------------------------*/
/*			MOVE ode function declations to myNP.h									 		*/
/*------------------------------------------------------------------------------------------*/

// 1st order ODE




// 2nd order ODE








int main(int argc, char* argv[])
{
	/*------------------------------------------------------------------------------------------*/
		/*==========================================================================*/
		/*					Part 1:  1st order ODE-IVP (single equation)			*/
		/*--------------------------------------------------------------------------*/
		/*   - A simple RC circuit response to sinusoidal voltage input				*/
		/*   - You can change the variable names									*/
		/*   - You can create additional functions									*/
		/*==========================================================================*/


		/*==========================================================================*/
		/*						Initialization										*/
		/*==========================================================================*/
		// Initial Conditions
	double a = 0;
	double b = 0.1;
	double h = 0.001;


	unsigned int N = ((b - a) / h) + 1;
	double y_EU[200] = { 0 };				//Cannot use sytanx of  y_EU[N]={0};	
	double y_RK2[200] = { 0 };
	double y_RK3[200] = { 0 };

	// IVP initial values
	double y_init = 0;


	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/
	// Euler's Explicit Method (given)
	odeEU(y_EU, odeFunc_rc, a, b, h, y_init); 

	//  Runge-Kutta 2nd order
	odeRK2(y_RK2, odeFunc_rc, a, b, h, y_init);

	//  Runge-Kutta 3rd order
	odeRK3(y_RK3, odeFunc_rc, a, b, h, y_init);



	/*==========================================================================*/
	/*							  Print your results							*/
	/*==========================================================================*/

	// Print outputs
	printf("----------------------------------------------------------------\n");
	printf("			       1st ODE - IVP  Results						\n");
	printf("----------------------------------------------------------------\n");

	printf("i\t t\t\t yEU\t\t yRK2\t\t yRK3\t \n\n");
	for (int i = 0; i < N; i++)
		printf("%d\t %0.5f\t %.5f\t %.5f\t %.5f\t \n", i, a + i * h, y_EU[i], y_RK2[i], y_RK3[i]);
	printf("\n");


	/*------------------------------------------------------------------------------------------*/

		/*==========================================================================*/
		/*					Part 2:  2nd order ODE-IVP								*/
		/*--------------------------------------------------------------------------*/
		/*   - A simple m-c-k  response to sinusoidal force input					*/
		/*   - You can change the variable names									*/
		/*   - You can create additional functions									*/
		/*==========================================================================*/

		/*==========================================================================*/
		/*						Initialization										*/
		/*==========================================================================*/
		//Parameter Definitions
	double t0 = 0.0 ;
	double tf = 1.0 ;
	h = 0.01 ;
	N = (tf - t0) / h + 1 ;

	double Y[200] = { 0 } ;
	double V[200] = { 0 } ;

	double Y2[200] = { 0 };
	double V2[200] = { 0 };

	// Initial values
	double y_initQ2 = 0 ;
	double v_initQ2 = 0.2 ;


	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/
	//  Runge-Kutta 2nd order
	sys2RK2(Y, V, odeFunc_mck, t0, tf, h, y_initQ2, v_initQ2);


	// Extra: Runge-Kutta 4th order
	sys2RK4(Y2, V2, odeFunc_mck, t0, tf, h, y_initQ2, v_initQ2); 


	/*==========================================================================*/
	/*							  Print your results							*/
	/*==========================================================================*/
	//
	printf("----------------------------------------------------------------\n");
	printf("			       2nd order ODE - IVP  Results					\n");
	printf("----------------------------------------------------------------\n");

	printf("i\t t\t\t Y\t\t V\t\t \n\n");

	for (int i = 0; i < N; i++)
		printf("%d\t %0.5f\t %.5f\t %.5f\t \n", i, a + i * h, Y[i], V[i]);


	printf("----------------------------------------------------------------\n");
	printf("			       4nd order ODE (EXTRA)- IVP  Results					\n");
	printf("----------------------------------------------------------------\n");


	printf("i\t t\t\t Y\t\t V\t\t \n\n");

	for (int i = 0; i < N; i++)
		printf("%d\t %0.5f\t %.5f\t %.5f\t \n", i, a + i * h, Y2[i], V2[i]); 


	FILE* fp1;
	fopen_s(&fp1, "yEU.txt", "w");

	for (int i = 0; i < 100; i++) {
		fprintf(fp1, "%f\n", y_EU[i]);
	}

	fclose(fp1);



	FILE* fp2;
	fopen_s(&fp2, "yRK2.txt", "w"); 

	for (int i = 0; i < 100; i++) {
		fprintf(fp2, "%f\n", y_RK2[i]);
	}

	fclose(fp2) ;
	

	FILE* fp3;
	fopen_s(&fp3, "yRK3.txt", "w");

	for (int i = 0; i < 100; i++) {
		fprintf(fp3, "%f\n", y_RK3[i]);
	}

	fclose(fp3);


	FILE* fp4;
	fopen_s(&fp4, "2RK2_Y.txt", "w");

	for (int i = 0; i < 100; i++) {
		fprintf(fp4, "%f\n", Y[i]);
	}

	fclose(fp4);


	FILE* fp5;

	fopen_s(&fp5, "2RK2_V.txt", "w");

	for (int i = 0; i < 100; i++) {
		fprintf(fp5, "%f\n", V[i]);
	}

	fclose(fp5);
	

	FILE* fp6 ;

	fopen_s(&fp6, "2RK4_Y.txt", "w") ;

	for (int i = 0; i < 100; i++) {
		fprintf(fp6, "%f\n", Y2[i]) ;
	}

	fclose(fp6) ;



	FILE* fp7;

	fopen_s(&fp7, "2RK4_V.txt", "w");

	for (int i = 0; i < 100; i++) {
		fprintf(fp7, "%f\n", V2[i]);
	}

	fclose(fp7);


	system("pause");
	return 0;
}


/*------------------------------------------------------------------------------------------*/



/*==========================================================================*/
/*							  Function Definitions							*/
/*==========================================================================*/

// Gradient function for ODE - 1st order 
double odeFunc_rc(const double t, const double v)
{
	// Input:	 y, t
	// Output:	 dydt 

	// system modeling parameters 
	double tau = 0.01;
	double f = 100;
	double Vm = 1;
	double omega = 2 * PI * f;
	double dvdt = 0;

	dvdt = - v / tau + Vm / tau * cos(omega*t) ; 

	return dvdt ;
}


// Gradient function for ODE - 2nd order 
void odeFunc_mck(double dYdt[], const double t, const double Y[])
{
	// Input:	 vecY = [y; v], t 
	// Output:	 dYdt = [dydt; dvdt] 

	// system modeling parameters 
	// [YOUR CODE GOES HERE]	
	// [YOUR CODE GOES HERE]	

	double m = 1.0 ; 
	double k = 6.9 ; 
	double c = 10 * 0.7 ;

	double FinDC = 2.0 ;    
	double f = 5
		.0 ;
	double Fin = FinDC * cos(2 * PI * f * t) ;



	// output
	// [YOUR CODE GOES HERE]	
	// dYdt[0] =_______________

	dYdt[0] = Y[1];

	// [YOUR CODE GOES HERE]	
	// dYdt[1] =_______________

	dYdt[1] = 1 / m * (Fin - c * Y[1] - k * Y[0]);

}





/*------------------------------------------------------------------------------------------*/
/*			MOVEfunction definitions to  myNP.cpp								 		*/
/*------------------------------------------------------------------------------------------*/





