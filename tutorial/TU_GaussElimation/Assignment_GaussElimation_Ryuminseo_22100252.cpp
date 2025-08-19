/*-------------------------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 19-10-2023
Language/ver     : C++ in MSVS2019

Description      : Tutorial code for using matrix structure
-------------------------------------------------------------------------------*/

#define ASGN		5		// enter your assignment number
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


	// Q1. 

	int rows_1 = matrix_A_1.rows;
	int cols_1 = matrix_A_1.cols;


	Matrix mat_answer_1 = createMat(rows_1, cols_1);
	Matrix vec_answer_1 = createMat(cols_1,1);


	gaussElim_Basic(matrix_A_1, vector_b_1, mat_answer_1, vec_answer_1);

	Matrix x_answer_1 = createMat(cols_1, 1);

	backsub(mat_answer_1, vec_answer_1, x_answer_1);

	printMat(matrix_A_1, "matA_Q1");
	printMat(vector_b_1, "vecb_Q1 ");
	printMat(mat_answer_1, "matU_Q1") ; 
	printMat(vec_answer_1, "vecd_Q1 ") ; 
	printMat(x_answer_1,   "vecx_Q1") ; 



	// Q2. 


	int rows_2 = matrix_A_2.rows;
	int cols_2 = matrix_A_2.cols;


	Matrix mat_answer_2 = createMat(rows_2, cols_2);
	Matrix vec_answer_2 = createMat(cols_2, 1);


	gaussElim_Basic(matrix_A_2, vector_b_2, mat_answer_2, vec_answer_2);

	Matrix x_answer_2 = createMat(cols_2, 1);

	backsub(mat_answer_2, vec_answer_2, x_answer_2);

	printMat(matrix_A_2, "matA_Q2");
	printMat(vector_b_2, "vecb_Q2 ");
	printMat(mat_answer_2, "matU_Q2");
	printMat(vec_answer_2, "vecd_Q2");
	printMat(x_answer_2, "vecx_Q2");



	// Q3.



	
	Matrix matrix_A_extra = txt2Mat(path, "mat_extra");

	int rows_extra = matrix_A_extra.rows;
	int cols_extra = matrix_A_extra.cols;

	Matrix L = eye(rows_extra, cols_extra);
	Matrix P = eye(rows_extra, cols_extra);
	Matrix U = createMat(rows_extra, cols_extra);

	Matrix vect_b = txt2Mat(path, "vect_extra");
	Matrix vect_x = createMat(rows_extra, 1);

	/*for (int i = 0; i < rows_extra; i++) {
		vect_b.at[i][1] = i + 1;
	}*/

	//gaussElim(matrix_A_extra, vector_b_extra, mat_answer_extra, vec_answer_extra, P);
	LUdecomp(matrix_A_extra, U, L, P);

	printMat(U, "U"); 
	printMat(P, "P"); 
	printMat(L, "L"); 

	Matrix vect_y = multMat(U, vect_x);

	fwdsub(L, multMat(P, vect_b),vect_y) ;
	backsub(U, vect_y, vect_x) ; 

	printMat(matrix_A_extra, "matA_extra") ; 
	printMat(vect_b, "b") ; 
	printMat(U, "U") ; 
	printMat(P, "P") ;  
	printMat(L, "L") ;  

	printMat(vect_y, "Y");
	printMat(vect_x, "x") ;  
	printMat(multMat(P, vect_b), "test") ; 
	




	// LU




	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/
	freeMat(matrix_A_1) ;		 
	freeMat(vector_b_1) ;      
	freeMat(matrix_A_2) ;      
	freeMat(vector_b_2) ;   


	freeMat(mat_answer_1) ;  
	freeMat(vec_answer_1) ;  
	freeMat(x_answer_1) ;    

	freeMat(mat_answer_2);
	freeMat(vec_answer_2);
	freeMat(x_answer_2);

	freeMat(P);
	freeMat(matrix_A_extra);
	freeMat(U);
	freeMat(L);

	freeMat(vect_b);
	freeMat(vect_x);
	freeMat(vect_y);


	// free other  created  matrices
	

	// test 

	/*Matrix ex1 = eye(5, 5);
	Matrix ex2 = txt2Mat(path, "prob1_matA");
	printMat(ex1, "test ex1");
	printMat(ex2, "test ex2");
	printMat(multMat(ex1,ex2),"test");
	freeMat(ex1);
	freeMat(ex2);*/

	system("pause");
	return 0;
}