/*------------------------------------------------------------------------------------------\
@ Numerical Programming  by Young-Keun Kim - Handong Global University

Author          : Young-Keun Kim
Created         : 01-04-2019
Modified        : 01-04-2024
Language/ver	: C in MSVS2017
Course			: Numerical Programming

Description     : Assignment 9 Curvefitting
/------------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------*/
/*			MODIFY Path and name of Headerfiles										 		*/
/*------------------------------------------------------------------------------------------*/
#include "../../include/myNP_22100252.h" 
#include "../../include/myMatrix_22100252.h" 




/*------------------------------------------------------------------------------------------*/
/*			MOVE the followings  to  myNP.h									 		*/
/*------------------------------------------------------------------------------------------*/
// Calculates coefficients of least squares regression - Line

// Calculates coefficients of least squares regression - Nth order polynomial


// Option 2:. Delete the below if you selected to use 1D arrays




//void expFit(double vecZ[], double vecX[], double vecY[]);
// or
//Matrix	expFit_mat(Matrix _X, Matrix _Y);


/*------------------------------------------------------------------------------------------*/
/*			MOVE the following  to  myMatrix.h									 		*/
/*------------------------------------------------------------------------------------------*/
// Create a Matrix from 1D-array
  // <-- move this to myMatrix.h, myMatrix.cpp 


int main(int argc, char* argv[])
{
	/*------------------------------------------------------------------------------------------*/
		/*==========================================================================*/
		/*					Part 1-1:  Polyfit(line)								*/
		/*--------------------------------------------------------------------------*/
		/*   - order n=1 linearFit													*/
		/*==========================================================================*/


		/*==========================================================================*/
		/*						Initialization										*/
		/*==========================================================================*/

		// Initial Conditions
	double T[] = { 30, 40, 50, 60, 70, 80 };
	double P[] = { 1.05, 1.07, 1.09, 1.14, 1.17, 1.21 };
	double Z_Q1[2] = { 0 } ; 
	int n = 1;	// nth order 
	int m_Q1 = 6;	// length of dataset 



	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/
	// Option 1: using 1D array
	n = 1 ;
	// [YOUR CODE GOES HERE]
	// [YOUR CODE GOES HERE]

	// Option 2	
	// Delete the below if you selected Option 1
	Matrix matT = arr2Mat(T, m_Q1, 1);
	Matrix matP = arr2Mat(P, m_Q1, 1);
	Matrix vecZ_Q1 = polyFit_mat(matT, matP, n);
	printMat(vecZ_Q1, "Z_Q1");



	/*------------------------------------------------------------------------------------------*/
		/*==========================================================================*/
		/*					Part 1-2:  Polyfit(nth )								*/
		/*--------------------------------------------------------------------------*/
		/*   - order n=2 or higher, polyfit											*/
		/*==========================================================================*/


		/*==========================================================================*/
		/*						Initialization										*/
		/*==========================================================================*/

		// Initial Conditions
	int m_Q2 = 16;				// data length
	
	double Stress[] = { 0, 3, 4.5, 5.8, 5.9, 5.8, 6.2, 7.4, 9.6, 15.6, 20.7, 26.7,31.1, 35.6, 39.3, 41.5 };
	
	double Strain[16] = { 0 } ;

	for (int k = 0; k < m_Q2; k++) {
		Strain[k] = 0.4 * k;
	}
	double Z_Q2[5] = { 0 };		// 4th order - a0 to a4


	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/
	// Option 1: using 1D array
	n = 4;	// nth order
	printf("Q2 \n");


	// Print results
	// [YOUR CODE GOES HERE]
	// [YOUR CODE GOES HERE]



	// Option 2	
	// Delete the below if you selected Option 1
	Matrix matStrain = arr2Mat(Strain, m_Q2, 1);
	Matrix matStress = arr2Mat(Stress, m_Q2, 1);
	//

	Matrix vecZ_Q2 = polyFit_mat(matStrain, matStress, n);
	printMat(vecZ_Q2, "Z_Q2");

	Matrix matStrain2 = arr2Mat(Strain, m_Q2, 1);
	Matrix matStress2 = arr2Mat(Stress, m_Q2, 1);
	Matrix vecZ_Q2_count = polyFit_mat_count(matStrain2, matStress2, n);
	printMat(vecZ_Q2_count, "Z_Q2 count");




	/*------------------------------------------------------------------------------------------*/
		/*==========================================================================*/
		/*					Part 2:  Exponential fit								*/
		/*--------------------------------------------------------------------------*/
		/*   yhat=a0*exp(a1*t)														*/
		/*   														*/
		/*==========================================================================*/

		/*==========================================================================*/
		/*						Initialization										*/
		/*==========================================================================*/
		// [YOUR CODE GOES HERE]
		// [YOUR CODE GOES HERE]
		// [YOUR CODE GOES HERE]

	n = 4;
	int m_Q3 = 15;
	double X_value[15] = { 0 }; 

	for (int i = 0; i < 15; i++) {
		X_value[i] = (i + 1) * 2;
	}

	double Y_value[15] = { 9.7, 8.1, 6.6, 5.1, 4.4, 3.7, 2.8, 2.4, 2.0, 1.6, 1.4, 1.1, 0.85, 0.69, 0.6 };
	double Y_value_log[15] = { 0 };

	for (int i = 0; i < m_Q3; i++) {
		Y_value_log[i] = log(Y_value[i]);
	}

	/*for (int i = 0; i < m_Q3; i++) {
		printf("Y_value_log : %f \n", Y_value_log[i]);
	}*/

	Matrix X_value_mat = arr2Mat(X_value, m_Q3, 1) ; 
	Matrix Y_value_mat = arr2Mat(Y_value_log, m_Q3, 1) ;
	Matrix vecZ_Q3 = expFit(X_value_mat, Y_value_mat) ; 

	//printf("C : %f \n", -1 * (1.0 / (vecZ_Q3.at[1][0] * 5*pow(10, 6))));
	//printMat(vecZ_Q3, "[ Z_Q3 ]") ;  


	system("pause");
	return 0;
}




