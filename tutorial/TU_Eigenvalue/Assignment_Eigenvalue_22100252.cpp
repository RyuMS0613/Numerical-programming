/*------------------------------------------------------------------------------------------\
@ Numerical Programming  by Young-Keun Kim - Handong Global University

Author          : Young-Keun Kim
Created         : 01-04-2019
Modified        : 01-04-2023
Language/ver	: C in MSVS2017
Course			: Numerical Programming

Description     : Assignment 7 Eigenvalue,Eigenvector
/------------------------------------------------------------------------------------------*/

#define ASGN	7		// enter your assignment number
#define EVAL	0		// set 0

// [ MODIFY Path and name of Headerfiles]
#include "../../include/myMatrix_22100252.h"
#include "../../include/myNP_22100252.h" 


// [ CREATE FUNCTIONS ]
// These must be in myMatrix.h, myMatrix.cpp
// Matrix eigval(Matrix A);
// void QRdecomp(Matrix A, Matrix Q, Matrix R);
// Matrix eigvec(Matrix A);
// void eig(Matrix A, Matrix V, Matrix D);




int main(int argc, char* argv[])
{
#if _WIN64 | _WIN32
	/*	 [воик DO NOT EDIT !!!]   Resources file path setting for evaluation	*/
	std::string path = "../../NP_Data/Assignment" + std::to_string(ASGN) + "/";
#elif __APPLE__
	//std::string path = "~/NP_Data/Assignment" + std::to_string(ASGN) + "/";
	std::string path = "../../NP_Data/Assignment" + std::to_string(ASGN) + "/";
#endif
#if EVAL
	path += "eval/";
#endif

	/*==========================================================================*/
	/*					Variables declaration & initialization					*/
	/*--------------------------------------------------------------------------*/
	/*   - You can change the variable names									*/
	/*   - However, you must use the specified file name						*/
	/*		   : For each assignment, the file name will be notified			*/
	/*==========================================================================*/
	Matrix matA_Q1 = txt2Mat(path, "prob1_matA");
	Matrix eigVals_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);
	Matrix eigVecs_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);
	Matrix matD_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);
	Matrix matV_Q1 = createMat(matA_Q1.rows, matA_Q1.cols);

	//printMat(matV_Q1, "[ matV_Q1 ]");
	//printMat(matD_Q1, "[ matD_Q1 ]");


	Matrix matA_Q2 = txt2Mat(path, "prob2_matA") ;
	Matrix eigVals_Q2 = createMat(matA_Q2.rows, matA_Q2.cols) ; 
	Matrix eigVecs_Q2 = createMat(matA_Q2.rows, matA_Q2.cols) ;
	Matrix matD_Q2 = createMat(matA_Q2.rows, matA_Q2.cols) ;
	Matrix matV_Q2 = createMat(matA_Q2.rows, matA_Q2.cols) ;



	//printMat(matV_Q2, "[ matV_Q2 ]");
	//printMat(matD_Q2, "[ matD_Q2 ]");




	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/

	// [Eigenval/vectors for Q1]
	//eigVals_Q1 = eigval(matA_Q1);
	//eigVecs_Q1 = eigvec(matA_Q1);

	//eigVals_Q2 = eigval(matA_Q2);
	//eigVecs_Q2 = eigvec(matA_Q2);

	//eig(matA_Q1, matV_Q1, matD_Q1);


	// [Eigenval/vectors for Q2]
	// YOUR CODE GOES HERE
	// YOUR CODE GOES HERE


	/*==========================================================================*/
	/*							  Print your results						*/
	/*==========================================================================*/

	// [Eigenval/vectors for Q1]
	printf("----------------------------------------------------------------------------------------------\n");
	printf("			       Eigenvalue & Eigenvector Results					\n");
	printf("----------------------------------------------------------------------------------------------\n");

	printf("\n----------------Q1----------------\n") ; 


	eig(matA_Q1, matV_Q1, matD_Q1) ; 

	printf("\n----------------Q2----------------\n") ; 

	eig(matA_Q2, matV_Q2, matD_Q2) ; 


	printf("\n[Eigen value of Q1]\n\n") ;
	printMat(eigVals_Q1, "[ Values_Q1 ]") ;
	printf("\n[Eigen vector]\n\n") ;
	printMat(eigVecs_Q1, "[ Vectors_Q1 ]") ;
	printf("\n\n") ;

	//// [Eigenval/vectors for Q2
	//printf("\n[Eigen value of Q2]\n\n");
	//printMat(eigVals_Q2, "[ Values_Q2 ]");
	//printf("\n[Eigen vector]\n\n");
	//printMat(eigVecs_Q2, "[ Vectors_Q2 ]");
	//printf("\n\n");


	//QRdecomp(matA_Q2,Q_2,R_2); 
	/*eigVals_Q2 = eigval(matA_Q2);
	eigVecs_Q2 = eigvec(matA_Q2);*/






	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/

	// [Deallocate for Q1]
	//freeMat(matA_Q1);
	////freeMat(valEig_Q1);	freeMat(vecEig_Q1);
	//freeMat(matD_Q1); 
	//freeMat(matV_Q1);

	//freeMat(matA_Q1); 
	//freeMat(matD_Q1);
	//freeMat(matV_Q1); 

	//freeMat(matA_Q2); 
	//freeMat(matD_Q2);
	//freeMat(matV_Q2);






	// [Deallocate for Q2]
	// YOUR CODE GOES HERE
	// YOUR CODE GOES HERE


	system("pause");
	return 0;
}