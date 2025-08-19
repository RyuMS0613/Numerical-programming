/*-------------------------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 19-10-2023
Language/ver     : C++ in MSVS2019

Description      : Tutorial code for using matrix structure
-------------------------------------------------------------------------------*/

#define ASGN		6		// enter your assignment number
#define EVAL		0		// [воик DO NOT EDIT !!!]


#include "../../include/myMatrix_22100252.h"

#include "../../include/myNP_22100252.h"


int main(int argc, char* argv[])
{

	/*	 [воик DO NOT EDIT !!!]   Resources file path setting for evaluation	*/
// #if _WIN64 | _WIN32
// 	std::string path = "C:/NP_Data/Assignment" + std::to_string(ASGN) + "/";

// #elif __APPLE__
// 	std::string path = "~/NP_Data/Assignment" + std::to_string(ASGN) + "/";
// #endif

	std::string path = "../../NP_Data/Assignment" + std::to_string(ASGN) + "/";
#if EVAL
	path += "eval/";
#endif



	// Read from datafile
	Matrix matrix_A_1 = txt2Mat(path, "prob1_matA");
	Matrix vector_b_1 = txt2Mat(path, "prob1_vecb");

	Matrix matrix_A_2 = txt2Mat(path, "prob2_matK");
	Matrix vector_b_2 = txt2Mat(path, "prob2_vecf");

	Matrix checking = txt2Mat(path, "prob3_mat");


	// checking

	int rows_0 = checking.rows;
	int cols_0 = checking.cols;


	Matrix P_0 = eye(rows_0, cols_0); 
	Matrix L_0 = eye(rows_0, cols_0); 
	Matrix U_0 = createMat(rows_0, cols_0);

	//scaled_LUdecomp(checking,U_0,L_0,P_0);
	//printMat(multMat(L_0, U_0),"LU");



	// Q1. 

	printf("----------------------------------------------\n");
	printf("--------------------- Q1 ---------------------\n");
	printf("----------------------------------------------\n");
	printf("\n");

	int rows_1 = matrix_A_1.rows ; 
	int cols_1 = matrix_A_1.cols ; 


	Matrix P_1 = eye(rows_1, cols_1) ; 
	Matrix L_1 = eye(rows_1, cols_1) ; 
	Matrix U_1 = createMat(rows_1, cols_1) ; 

	Matrix x = createMat(rows_1, 1); 

	//scaled_LUdecomp(matrix_A_1,U_1,L_1,P_1) ; 

	

	// Q2. 

	printf("\n");
	printf("----------------------------------------------\n");
	printf("--------------------- Q2 ---------------------\n");
	printf("----------------------------------------------\n");
	printf("\n");

	int rows_2 = matrix_A_2.rows ; 
	int cols_2 = matrix_A_2.cols ; 
	 
	Matrix P_2 = eye(rows_2,cols_2) ; 
	Matrix L_2 = eye(rows_2, cols_2) ; 
	Matrix U_2 = createMat(rows_2, cols_2) ; 

	Matrix x2 = createMat(rows_2, 1);

	

	//scaled_LUdecomp(matrix_A_2,U_2,L_2,P_2) ; 




	// Problem 2.

	printf("\n");
	printf("-----------------------------------------------\n");
	printf("------------------ Problem 2 ------------------\n");
	printf("-----------------------------------------------\n");
	printf("\n");


	Matrix L_3 = eye(rows_1, cols_1);
	Matrix P_3 = eye(rows_1, cols_1);
	Matrix U_3 = createMat(rows_1, cols_1);

	//Q1 
	//LUdecomp(matrix_A_1, U_3, L_3, P_3); 
	//LUsolve(L_3,U_3,P_3,vector_b_1,x); 
	
	//Q2
	//LUdecomp(matrix_A_2, U_2, L_2, P_2);
	//LUsolve(L_2, U_2, P_2, vector_b_2, x2);


	printf("\n");
	printf("----------------------------------------------\n");
	printf("------------------ Problem 3 ------------------\n");
	printf("----------------------------------------------\n");
	printf("\n");



	Matrix Ainv1 = eye(rows_1, cols_1) ; 
	Matrix Ainv2 = eye(rows_2, cols_2) ; 

	// print invert reslut

	printMat(matrix_A_1,"Q1");
	invMat(matrix_A_1, Ainv1);
	printMat(multMat(Ainv1, vector_b_1), "xvalue_Q1");

	printMat(matrix_A_2, "Q2");
	invMat(matrix_A_2, Ainv2);
	printMat(multMat(Ainv2, vector_b_2), "xvalue_Q2");


	




	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/
	
	freeMat(matrix_A_1);
	freeMat(vector_b_1);
	freeMat(matrix_A_2);
	freeMat(vector_b_2);

	freeMat(checking);
	freeMat(P_0);
	freeMat(U_0);
	freeMat(L_0);


	// Q1 free 

	freeMat(P_1);
	freeMat(U_1);
	freeMat(L_1);
	freeMat(x);

	// Q2 free

	freeMat(P_2); 
	freeMat(U_2); 
	freeMat(L_2); 
	freeMat(x2);


	// Q3 free
	freeMat(P_3);
	freeMat(U_3);
	freeMat(L_3);


	freeMat(Ainv1);
	freeMat(Ainv2);


	system("pause");
	return 0;
}