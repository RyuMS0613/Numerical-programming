/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 14-10-2024
Language/ver     : C++ in MSVS2019

Description      : myMatrix.cpp
----------------------------------------------------------------*/

#include "myMatrix_22100252.h"



// Free a memory allocated matrix
void	freeMat(Matrix _A)
{
	// 1. Free allocated column memory
	for (int i = 0; i < _A.rows; i++) {
		free(_A.at[i]);
		}
	// 2. Free allocated row memory
	free(_A.at);
}

// Create a matrix from a text file
Matrix	txt2Mat(std::string _filePath, std::string _fileName)
{
	std::ifstream file;
	std::string temp_string, objFile = _filePath + _fileName + ".txt";
	int temp_int = 0, nRows = 0;

	file.open(objFile);
	if (!file.is_open()) {
		printf("\n*********************************************");
		printf("\n  Could not access file: 'txt2Mat' function");
		printf("\n*********************************************\n");
		return createMat(0, 0);
	}
	while (getline(file, temp_string, '\t'))
		temp_int++;
	file.close();

	file.open(objFile);
	while (getline(file, temp_string, '\n'))
		nRows++;
	file.close();

	int nCols = (temp_int - 1) / nRows + 1;
	Matrix Out = createMat(nRows, nCols);

	file.open(objFile);
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < nCols; j++) {
			file >> temp_string;
			Out.at[i][j] = stof(temp_string);
		}
	file.close();

	return Out;
}

// Create Matrix with specified size
Matrix	createMat(int _rows, int _cols)
{
	// check matrix dimension
	if (_rows < 0 || _cols < 0) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'createMat' function");
		printf("\n****************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out;
	// 1. Allocate row array first
	Out.at = (double**)malloc(sizeof(double*) * _rows);
	// 2. Then, allocate column 
	for (int i = 0; i < _rows; i++)
		Out.at[i] = (double*)malloc(sizeof(double) * _cols);
	// 3. Initialize row & column values of a matrix
	Out.rows = _rows;
	Out.cols = _cols;

	// 4. Initialize with zero (optional)
	initMat(Out, 0);
	return Out;
}


// initialization of Matrix elements
void	initMat(Matrix _A, double _val)
{
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			_A.at[i][j] = _val;
}

// Print matrix
void	printMat(Matrix _A, const char* _name)
{
	printf("%s =\n", _name);
	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++)
			printf("%15.4f\t", _A.at[i][j]);
		printf("\n");
	}
	printf("\n");
}

// Matrix addition
Matrix	addMat(Matrix _A, Matrix _B)
{
	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = createMat(_A.rows, _B.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j] = _A.at[i][j] + _B.at[i][j];

	return Out;
}






//////////////////////////////////////////////////////////////////
/*				Tutorial	&  Assignment						*/
//////////////////////////////////////////////////////////////////

// Create matrix of all zeros
extern Matrix zeros(int _rows, int _cols)
{
	Matrix Out = createMat(_rows, _cols);
	initMat(Out, 0);
	return Out;
}


// Create matrix of all ones
extern Matrix	ones(int _rows, int _cols)
{
	Matrix Out = createMat(_rows, _cols);
	initMat(Out,1.0);

	return Out;
}


// Create identity matrix
extern Matrix eye(int _rows, int _cols)
{
	Matrix Out = createMat(_rows, _cols);

	for (int i = 0; i < _rows; i++) {
		for (int j = 0; j < _cols; j++) {
			if (i == j) Out.at[i][j] = 1.0;
		
		}
	}

	return Out;
}


// Matrix subtraction
extern	Matrix	subMat(Matrix _A, Matrix _B) {
	Matrix Out = createMat(_A.rows, _A.cols);
	
	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++) {

			 Out.at[i][j] = _A.at[i][j]- _B.at[i][j];

		}
	}

	return Out;
}


// Multiply  matrix A and matrix B  OUT=AB
extern	Matrix	multMat(Matrix _A, Matrix _B) 
{
	Matrix Out = createMat(_A.rows, _B.cols);

	if (_A.cols != _B.rows) {
		printf("Error");
	}
	else {
	
		for (int i = 0; i < Out.rows; i++) {
			for (int j = 0; j < Out.cols; j++) {

				for (int k = 0; k < _B.rows; k++) {

					Out.at[i][j] = Out.at[i][j] + _A.at[i][k] * _B.at[k][j];
				}

			}
		}
	
	}

	return Out;

}


extern	void multMat_void(Matrix _A, Matrix _B, Matrix _AB) {

	initMat(_AB,0.0); 

	if (_A.cols != _B.rows) {
		printf("Error");
	}
	else {

		for (int i = 0; i < _A.rows; i++) {

			for (int j = 0; j < _B.cols; j++) {

				for (int k = 0; k < _B.rows; k++) {

					_AB.at[i][j] = _AB.at[i][j] + _A.at[i][k] * _B.at[k][j] ; 

				}

			}
		}

	}
}


// Multiply  matrix A with a scalar k
extern	Matrix	smultMat(Matrix _A, double _k) {
	Matrix Out = createMat(_A.rows, _A.cols);

	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++) {

			Out.at[i][j] = _k * _A.at[i][j] ;

		}
	}
	return Out;
}


// Create Transpose matrix
extern	Matrix	transpose(Matrix _A) {
	Matrix Out = createMat(_A.cols, _A.rows);

	for (int i = 0; i < _A.cols; i++) {
		for (int j = 0; j < _A.rows; j++) {

			Out.at[i][j] =  _A.at[j][i];

		}
	}	

	return Out;
}


// Copy matrix
void copyMat(Matrix _A, Matrix _B) {

	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++) {

			_A.at[i][j] = _B.at[i][j] ;

		}
	}
}



// gaussElimination
void gaussElim_Basic(Matrix A, Matrix b, Matrix U, Matrix d) {

	double mult = 0.0;
	int m = A.rows;
	int n = A.cols;

	copyMat(U, A);
	copyMat(d, b);

	// Error check (정방행렬 아닌 경우)
	if (m != n) {
		initMat(U,0.0); // U를 0으로 초기화
		printf("정방행렬이 아닙니다");
	}

	else
	{
		for (int k = 0; k < m - 1; k++) {

			// Error check (Pivot이 0인 경우), 
			if (U.at[k][k] == 0) {
				printf("Pivot is 0 (zero)");
				initMat(U, 0.0); // 이전 계산내용 0으로 초기화
				break;
			}

			for (int i = k + 1; i < m; i++) {

				mult = U.at[i][k] / U.at[k][k];
				d.at[i][0] = d.at[i][0] - mult * d.at[k][0];

				for (int j = k; j < n; j++) {
					U.at[i][j] = U.at[i][j] - mult * U.at[k][j];
				}
			}
		}
	}
}



// gaussElimination with povoting
void gaussElim(Matrix A, Matrix b, Matrix U, Matrix d, Matrix P) {
	
	double mult = 0.0 ;
	double temp = 0.0 ; 

	int m = A.rows ;
	int n = A.cols ;

	double max = 0.0 ;
	int max_number = 0 ;

	copyMat(U, A); 
	copyMat(d, b);

	// Error check (정방행렬 아닌 경우)
	if (m != n) 
	{
		printf("We will only consider square matrix A (n by n) for this assignment");
		initMat(U, 0.0);
	}
	else 
	{
	
		for (int k = 0; k < m - 1; k++) {

			// Error check (Pivot이 0인 경우), 
			if (U.at[k][k] == 0) {
				printf("Pivot is 0 (zero)") ; 
				// gaussElim 함수를 그만 두기 전 zero로 바꾸어주고 난 후 끝냄. (이전 계산 내용 초기화)
				initMat(U, 0.0) ;
				break; 
			}


			// 가장 큰 pivot 구하기

			max = fabs(U.at[k][k]);

			for (int mn = k + 1 ; mn < m  ; mn++) {
				
				if (fabs(U.at[mn][k]) > max) {
					max = fabs(U.at[mn][k]);
					max_number = mn ; 
				}

			}

			// 순서 바꾸기 

			for (int i = 0 ; i < n ; i++) {
				temp = U.at[max_number][i];
				U.at[max_number][i] = U.at[k][i];
				U.at[k][i] = temp;

				temp = P.at[max_number][i];
				P.at[max_number][i] = P.at[k][i];
				P.at[k][i] = temp;
				
			}

			// Gauss Elimination

			for (int i = k + 1; i < m; i++) {

				mult = U.at[i][k] / U.at[k][k];

				d.at[i][0] = d.at[i][0] - mult * d.at[k][0];

				for (int j = k; j < n; j++) {

					U.at[i][j] = U.at[i][j] - mult * U.at[k][j];

				}

			}

		}
	}


}


void LUdecomp(Matrix A, Matrix U, Matrix L) {

	double mult = 0.0;
	double temp = 0.0;
	int m = A.rows;
	int n = A.cols;
	double max = 0.0;
	int max_index = 0;

	int check = 0;

	copyMat(U, A);

	printMat(A, "[ A ] ");

	// Error check (정방행렬 아닌 경우)
	if (m != n)
	{
		printf("We will only consider square matrix A (n by n) for this assignment");
		initMat(U, 0.0);
	}
	else
	{
		for (int k = 0; k < m - 1; k++) {


			// Gauss Elimination + L 행렬 값 입력 
			for (int i = k + 1; i < m; i++) {

				mult = U.at[i][k] / U.at[k][k];
				L.at[i][k] = mult;

				for (int j = k; j < n; j++) {

					U.at[i][j] = U.at[i][j] - mult * U.at[k][j];
				}
			}

			printf("At k = %d \n\n", k);
			printMat(L, "[ L ]");
			printMat(U, "[ U ]");

		}
	}
}


// LUdecomp with pivoting
void LUdecomp_pivoting(Matrix A, Matrix U, Matrix L, Matrix P) {

	double mult = 0.0 ;
	double temp = 0.0 ;
	int m = A.rows ; 
	int n = A.cols ; 
	double max = 0.0 ; 
	int max_index = 0 ;

	int check = 0;

	copyMat(U, A);

	printMat(A, "[ A ] ");

	// Error check (정방행렬 아닌 경우)
	if (m != n)
	{
		printf("We will only consider square matrix A (n by n) for this assignment");
		initMat(U, 0.0);
	}
	else
	{
		for (int k = 0; k < m - 1; k++) {

			check = 0 ;

			max = fabs(U.at[k][k]); 
			for (int i = k + 1; i < m; i++) { 

				if (fabs(U.at[i][k]) > max) {
					max = fabs(U.at[i][k]);
					max_index = i;
					check = 1;
				}
			}

			if (check == 0) {
				max_index = k;
			}

			// 순서 바꾸기
			for (int i = 0; i < n; i++) {
				temp = U.at[max_index][i];
				U.at[max_index][i] = U.at[k][i];
				U.at[k][i] = temp;

				temp = P.at[max_index][i];
				P.at[max_index][i] = P.at[k][i];
				P.at[k][i] = temp;
			}

			// L 행렬 순서 바꾸기
			for (int i = 1 ; i <= k; i++) { 
				temp = L.at[max_index][k - i]; 
				L.at[max_index][k - i] = L.at[k][k - i]; 
				L.at[k][k - i] = temp; 
			} 

			// Gauss Elimination + L 행렬 값 입력 
			for (int i = k + 1; i < m; i++) { 

				mult = U.at[i][k] / U.at[k][k] ; 
				L.at[i][k] = mult ; 

				for (int j = k; j < n; j++) {

					U.at[i][j] = U.at[i][j] - mult * U.at[k][j] ;
				}
			}

			/*printf("At k = %d \n\n",k);
			printMat(P,"[ P ]") ;
			printMat(L,"[ L ]") ;
			printMat(U,"[ U ]") ;*/

		}
	}
}


// LUdecomp with scaled partial pivoting  
void scaled_LUdecomp(Matrix A, Matrix U, Matrix L, Matrix P) {

	double mult = 0.0;
	double temp = 0.0;
	int m = A.rows;
	int n = A.cols;

	double max_cols = 0.0;
	double max_rows = 0.0;
	int max_index = 0;


	copyMat(U, A);

	printMat(U, "[ U ] ");

	// Error check (정방행렬 아닌 경우)
	if (m != n)
	{
		printf("We will only consider square matrix A (n by n) for this assignment");
		initMat(U, 0.0);
	}
	else
	{
		for (int k = 0; k < m-1 ; k++) {

			printf("\niteration: %d", k);
			printMat(P, "P");
	
			max_rows = 0;

			// scaled pivot 구하기 

			for (int i = k ; i < m ; i++) {

				max_cols = fabs(U.at[i][k]);
				for (int j = k ; j < n ; j++) {
					if (max_cols < fabs(U.at[i][j])) max_cols = fabs(U.at[i][j]); 					
				}

				if (max_cols == 0) {
					printf("0으로 나누어지면 안됩니다");
					break;
				}

				if (max_rows < fabs(U.at[i][k] / max_cols)) { 
					max_rows = fabs(U.at[i][k] / max_cols) ;
					max_index = i ; 
				} 
			
			} 


			// 순서 바꾸기
			for (int i = 0; i < n; i++) {
				temp = U.at[max_index][i];
				U.at[max_index][i] = U.at[k][i];
				U.at[k][i] = temp;

				temp = P.at[max_index][i];
				P.at[max_index][i] = P.at[k][i];
				P.at[k][i] = temp;
			}

			// L 행렬 순서 바꾸기
			for (int i = 1; i <= k; i++) {
				temp = L.at[max_index][k - i];
				L.at[max_index][k - i] = L.at[k][k - i];
				L.at[k][k - i] = temp;
			}



			// Gauss Elimination + L 행렬 값 입력
			for (int i = k + 1; i < m; i++) {

				mult = U.at[i][k] / U.at[k][k];
				L.at[i][k] = mult;

				for (int j = k; j < n; j++) {

					U.at[i][j] = U.at[i][j] - mult * U.at[k][j];
				}
			}

			printf("At k = %d \n\n", k);


			printMat(P, "[ P ]");
			printMat(L, "[ L ]");
			printMat(U, "[ U ]");

		}
	}
}



// solve Ax = b 
Matrix solve_system(Matrix A, Matrix b) {

	int n = A.rows; 
	
	Matrix U = createMat(n,n);
	Matrix L = eye(n,n);
	Matrix P = eye(n,n);

	Matrix x = zeros(n,1); 
	Matrix y = zeros(n,1);

	LUdecomp(A, U, L);
	fwdsub(L, b, y); 
	backsub(U, y, x);


	return x ; 

}

// solve Ax = b with pivoting
Matrix solve_system_pivoting(Matrix A, Matrix b) {

	int n = A.rows;

	Matrix U = createMat(n, n);
	Matrix L = eye(n, n);
	Matrix P = eye(n, n);

	Matrix x = zeros(n, 1);
	Matrix y = zeros(n, 1);

	LUdecomp_pivoting(A, U, L, P); 

	fwdsub(L, multMat(P, b), y); 
	backsub(U, y, x);

	printMat(U , "U") ;
	printMat(L , "L") ;
	printMat(y , "y") ;
	printMat(P , "P") ;

	return x;

}


// LU solve
void LUsolve(Matrix L, Matrix U, Matrix P, Matrix b ,Matrix x) {

	Matrix y = multMat(U, x);

	fwdsub(L, multMat(P, b), y);

	printMat(y, "[ y ]");
	backsub(U, y, x);

	printMat(x,"[ x ]");

	freeMat(y);

}


// backward sub
void backsub(Matrix U, Matrix d, Matrix x) {

	double sum = 0.0 ;

	int m = U.rows;
	int n = U.cols;

	if (U.at[0][0] != 0) {
		for (int i = m - 1; i >= 0; i--) {

			sum = 0.0;

			for (int j = i + 1; j < n; j++) { 

				sum = sum + U.at[i][j] * x.at[j][0];
			}

			x.at[i][0] = (d.at[i][0] - sum) / U.at[i][i];

		}
	}
	else {
		printf(" Pivot is Zero ! \n");

	}
}


// forward sub 
void fwdsub(Matrix L, Matrix d, Matrix x) {

	double sum = 0.0;

	int m = L.rows ;
	int n = L.cols ;

	if (L.at[0][0] != 0) {
		for (int i = 0 ; i < n; i++) {

			sum = 0.0 ; 

			for (int j = 0 ; j < i; j++) {

				sum = sum + L.at[i][j] * x.at[j][0];
			}

			x.at[i][0] = (d.at[i][0] - sum) / L.at[i][i];

		}
	}
	else {
		printf(" Pivot is Zero ! \n");
	
	}
}


// invert Matrix 
void invMat(Matrix A, Matrix Ainv) {

	double mult = 0.0 ; 
	int m = A.rows ; 
	int n = A.cols ; 

	double temp_iden = 0 ;

	Matrix iden = createMat(m,n) ; 
	copyMat(iden,A) ; 


	for (int k = 0 ; k < m ; k++) { 

		temp_iden = iden.at[k][k] ; 


		// povot 1 만들기 + 해당 열 1 만들기
		for (int i = 0 ; i < n ; i++) { 
			iden.at[k][i] = iden.at[k][i] / temp_iden  ; 
			Ainv.at[k][i] = Ainv.at[k][i] / temp_iden  ; 
		} 


		for (int i = 0 ; i < m ; i++) { 
			
			if (i == k) continue ;

			// mult 정하기 
			mult =	iden.at[i][k] ;

			// 행 빼주기
			for (int j = 0 ; j < n; j++) {
				iden.at[i][j] = iden.at[i][j] - mult * iden.at[k][j] ;  
				Ainv.at[i][j] = Ainv.at[i][j] - mult * Ainv.at[k][j] ; 
			}
		}


	
	}

	//printMat(iden, "[ iden ]");
	//printMat(Ainv, "[ Ainv ]");

}


// Eigenvalue && Eigenvector


Matrix eigval(Matrix A) {

	Matrix Q = eye(A.rows, A.cols) ;  
	Matrix U = createMat(A.rows, A.cols);

	Matrix value = createMat(A.rows, 1) ; 

	Matrix R = createMat(A.rows, A.cols) ;

	copyMat(R, A); 
	copyMat(U, A); 

	if (A.cols != A.rows) {
		printf("정방 행렬이 아닙니다.");

		initMat(value, 0.0);
	}
	else
	{

		for (int i = 0; i < 100; i++) {
			QRdecomp(U, Q, R) ;
			//printMat(Q, "Q");
			//printMat(R, "R");
			multMat_void(R, Q, U) ;
			//printMat(U, "U");
		}

		for (int i = 0; i < A.rows; i++) {
			for (int j = 0; j < A.cols; j++) {
				if (i == j) {
					value.at[i][0] = U.at[i][j];
				}
			}
		}

		
		//printMat(U,"U");
		double max = 0.0 ;
		int max_index = 0 ;
		double temp = 0.0 ;

		// 순서 바꾸기 (sorting) 
		for (int i = 0; i < A.rows; i++) {
			max = value.at[i][0];
			max_index = i;
			for (int j = i; j < A.rows; j++) {
				if (max < value.at[j][0]) {
					max_index = j;
				}
			}
			temp = value.at[i][0];
			value.at[i][0] = value.at[max_index][0];
			value.at[max_index][0] = temp;

			//printMat(value, "value") ;  
		}

	}
	return value ;
}



void QRdecomp(Matrix A, Matrix Q, Matrix R) {

	int n = A.rows ;
	double sum = 0.0 ;
	double C_norm = 0.0 ;
	double V_in = 0.0 ;

	Matrix C = createMat(n,1) ;
	Matrix e = createMat(n,1) ;
	Matrix V = createMat(n,1) ;
	Matrix VT ; 
	Matrix H = createMat(n,n) ;

	Matrix Q_temp = eye(n, n) ; 
	Matrix R_temp = createMat(n, n);

	copyMat(R_temp,A);

	for (int i = 0; i < n - 1; i++) {
	
		initMat(C,0.0);

		for (int j = i ; j < n; j++) {
			C.at[j][0] = R_temp.at[j][i] ;
		}

		initMat(e, 0.0);

		e.at[i][0] = 1;
		
		if (C.at[i][0] < 0) {
			e.at[i][0] = -1;
		}

		sum = 0.0 ;

		for (int j = 0; j < n; j++) {
			sum = sum + pow(C.at[j][0],2); 
		}

		C_norm = sqrt(sum);


		initMat(V, 0.0);

		for (int j = 0; j < n; j++) {
			V.at[j][0] = C.at[j][0] + C_norm * e.at[j][0] ; 
		}

		//printMat(V, "[ V ]");


		VT = multMat(V,transpose(V));
		
		V_in = 0.0 ;

		for (int j = 0; j < n; j++) {
			V_in = V_in + pow(V.at[j][0], 2);
		}

		if (V_in == 0) {
			printf("0으로 나눌 수 없습니다.");
			break;
		}

		for (int j = 0; j < n; j++) {
			
			for (int k = 0; k < n; k++) {

				if (j == k) {
					H.at[j][k] = 1 - 2 * ( VT.at[j][k] / V_in ) ;
				}
				else {
					H.at[j][k] = -2 * ( VT.at[j][k] / V_in );
				}		
			}
		}

		//printMat(H,"[ H ]");

		multMat_void(Q_temp,H,Q) ; 
		multMat_void(H,R_temp,R) ; 
		
		
		copyMat(Q_temp, Q) ; 
		copyMat(R_temp, R) ;

		//printMat(Q, "[ Q ]");
		//printMat(R, "[ R ]");

	}
}


Matrix eigvec(Matrix A) {

	Matrix value = createMat(A.rows, 1); 
	value = eigval(A);

	Matrix eigenvect = ones(A.rows,A.cols);

	Matrix B = createMat(A.rows,A.cols);
	Matrix B_1 = createMat(A.rows-1, A.cols-1);
	Matrix B_b = createMat(A.rows-1, 1);

	Matrix value_temp = createMat(A.rows-1, 1);
	Matrix I = eye(A.rows,A.cols);
	
	Matrix norm_temp = eye(A.rows,1);

	int B_1_rows = 0;
	int B_1_cols = 0;
	int B_b_rows = 0;

	int value_temp_rows = 0;
	double norm_v = 0.0;
	
	if (A.rows == 2) {

		for (int i = 0; i < A.rows; i++) {

			B = subMat(A, smultMat(I, value.at[i][0]));

			for (int j = 0; j < A.cols; j++) {
				if (i == j) continue;
				
				if (B.at[j][j] == 0) {
					printf("decvided by 0 \n");
					break;
				}
				eigenvect.at[j][i] = -1 * B.at[j][i] / B.at[j][j] ; 
			}

			for (int j = 0; j < A.cols; j++) {
				norm_temp.at[j][0] = eigenvect.at[j][i];
			}

			norm_v = norm_2(norm_temp);

			if (norm_v == 0){
				printf("decvided by 0 \n");
				break;
			}

			for (int j = 0; j < A.cols; j++) {
				eigenvect.at[j][i] = eigenvect.at[j][i] / norm_v;
			}

		}
	
	}
	else if (A.rows == 3) {

		//printMat(eigenvect, "[ eigenvect - 1 ]");

		for (int i = 0; i < A.rows; i++) {

			B = subMat(A,smultMat(I,value.at[i][0]));	

			//printMat(B, "B");

			Matrix B_1_inv = eye(A.rows-1, A.cols-1);


			// B_1_cols

			B_1_cols = 0 ;

			for (int j = 0; j < A.rows ; j++) {
				if (i == j) continue;
				B_1_rows = 0 ;
				for (int k = 0; k < A.cols ; k++) {
					if (i == k) continue;

					B_1.at[B_1_rows][B_1_cols] = B.at[j][k] ;
					B_1_rows++;
					
				}
				B_1_cols++ ;
			}
 
			
			// B-b 만들기
			B_b_rows = 0;

			for (int j = 0; j < A.rows; j++) {
				if (i == j) continue;
				B_b.at[B_b_rows][0] = -B.at[j][i] ;
				B_b_rows  ++;
			}

			//printMat(B_1, "B_1");

			invMat(B_1, B_1_inv);

			//printMat(B_1_inv, "B_1_inv");
			value_temp = multMat(B_1_inv,B_b); 

			//printMat(value_temp, "value_temp"); 

			value_temp_rows = 0 ; 

			for (int j = 0; j < A.cols; j++) { 
				if (i == j) continue ; 
				eigenvect.at[j][i] = value_temp.at[value_temp_rows][0] ; 
				value_temp_rows++ ; 
			} 


			for (int j = 0; j < A.cols; j++) { 
				norm_temp.at[j][0] = eigenvect.at[j][i] ; 
			} 

			norm_v = norm_2(norm_temp);

			if (norm_v == 0) {
				printf("decvided by 0 \n");
				break;
			}

			for (int j = 0; j < A.cols; j++) {
				eigenvect.at[j][i] = eigenvect.at[j][i] / norm_v ;
			}


			//printMat(eigenvect, "result"); 

		}
	}
	else {
		printf("Input only 2X2, 3x3 Matrix");
	}

	return eigenvect ;
}



double norm_2(Matrix A) {

	double out = sqrt(multMat(transpose(A), A).at[0][0]);

	return out;
}


void eig(Matrix A, Matrix V, Matrix D) {


	Matrix value_vect = createMat(A.rows, 1); 
	Matrix temp = eigvec(A); 
	copyMat(V, temp); 

	value_vect = eigval(A); 

	for (int i = 0; i < A.rows; i++) {

		for (int j = 0; j < A.cols; j++) {
			if (i == j)
				D.at[j][j] = value_vect.at[i][0];
		}
	}

	printMat(V, "[ V ]"); 
	printMat(D, "[ D ]"); 

}



Matrix	linearFit_mat(Matrix _vecX, Matrix _vecY) {

	// Initialization	
	double Sx = 0;
	double Sxx = 0;
	double Sxy = 0;
	double Sy = 0;
	double a1 = 0; 
	double a0 = 0; 


	// Check m = length(X) and length(Y)
	int mx = _vecX.rows;
	int my = _vecY.rows;

	if (mx != my) {
		printf("ERROR: The number of elements in x must be the same as in y.");
	}
	else {

		for (int i = 0; i < mx; i++) {
			Sx = Sx + _vecX.at[i][0] ;
			Sy = Sy + _vecY.at[i][0] ;
			Sxx = Sxx + _vecX.at[i][0] * _vecX.at[i][0] ;
			Sxy = Sxy + _vecX.at[i][0] * _vecY.at[i][0] ;
		}

		a0 = (Sxx * Sy - Sx * Sxy) / (mx * Sxx - Sx * Sx) ; 
		a1 = (-Sx * Sy + Sxy * mx) / (mx * Sxx - Sx * Sx) ; 

	}

	/*printf("a0 : f% \n",a0) ;
	printf("a1 : f% \n",a1) ; */

	double z_array[2] = { a0, a1 };

	return arr2Mat(z_array, 2, 1); 
}


Matrix	arr2Mat(double* _1Darray, int _rows, int _cols)
{
	Matrix Output = createMat(_rows, _cols);

	for (int i = 0; i < _rows; i++)
		for (int j = 0; j < _cols; j++)
			Output.at[i][j] = _1Darray[i * _cols + j];

	return Output;
}


Matrix	polyFit_mat(Matrix _vecX, Matrix _vecY, int n) {

	// Initialization	
	Matrix vecZ = createMat(n + 1, 1);


	double a1 = 0.0 ;
	double a0 = 0.0 ;

	double Sx_temp = 0.0 ;
	double Sxy_temp = 0.0 ;

	Matrix Sx = createMat(2 * n + 1, 1); 
	Matrix Sxy = createMat(n + 1, 1); 
	Matrix S = createMat(n + 1, n + 1); 
	Matrix S_inv = eye(n + 1, n + 1); 

	int mx = _vecX.rows;
	int my = _vecY.rows;


	if (mx != my) {
		printf("ERROR: The number of elements in x must be the same as in y.");
	}

	if (n == 1)
	{
		vecZ= linearFit_mat(_vecX, _vecY) ; 

		return vecZ ;
	}
	else
	{

		for (int i = 0; i < 2 * n + 1; i++) {

			Sx_temp = 0;
			for (int j = 0; j < mx; j++) {
				Sx_temp = Sx_temp + pow(_vecX.at[j][0], i);
			}

			Sx.at[i][0] = Sx_temp ;
		}

		for (int i = 0; i < n + 1; i++) {
			for (int j = 0; j < n + 1; j++) {
				S.at[i][j] = Sx.at[i + j][0];
			}
		}


		for (int i = 0; i < n + 1; i++) {
			Sxy_temp = 0.0;
			for (int j = 0; j < mx; j++) {
				Sxy_temp = Sxy_temp + _vecY.at[j][0] * pow(_vecX.at[j][0], i);
			}
			Sxy.at[i][0] = Sxy_temp;
		}

		//printMat(S," S ");

		invMat(S, S_inv);

		multMat_void(S_inv, Sxy, vecZ);

	}

	
	return vecZ;
}


Matrix	polyFit_mat_count(Matrix _vecX, Matrix _vecY, int n) {

	int count = 0;

	// Initialization	
	Matrix vecZ = createMat(n + 1, 1);


	double a1 = 0.0;
	double a0 = 0.0;

	double Sx_temp = 0.0;
	double Sxy_temp = 0.0;

	Matrix Sx = createMat(2 * n + 1, 1);
	Matrix Sxy = createMat(n + 1, 1);


	Matrix S = createMat(n + 1, n + 1);
	Matrix St = createMat(n + 1, n + 1);
	Matrix D = createMat(n + 1, n + 1);

	Matrix result = createMat(n + 1, n + 1);


	Matrix S_inv = eye(n + 1, n + 1);


	int mx = _vecX.rows;
	int my = _vecY.rows;


	Matrix A = createMat(mx, n+1) ; 
	Matrix b = createMat(n+1, 1) ; 
	

	if (mx != my) {
		printf("ERROR: The number of elements in x must be the same as in y.");
	}


	if (n == 1)
	{
		vecZ = linearFit_mat(_vecX, _vecY) ;

		return vecZ ;
	}
	else
	{

		Sx.at[0][0] = mx ;
		D.at[0][0] = mx ;

		for (int i = 1; i < n+1; i++) {

			for (int j = 0; j < mx ; j++) {
				Sx.at[i * 2 - 1][0] = Sx.at[i * 2 - 1][0] + pow(_vecX.at[j][0], 2*i-1) ;
				Sx.at[i * 2][0] = Sx.at[i * 2][0] + pow(_vecX.at[j][0], 2*i) ;

				if (i == 1) { 
					b.at[0][0] = b.at[0][0] + _vecY.at[j][0];
				} 
				
				b.at[i][0] = b.at[i][0] + _vecY.at[j][0] * pow(_vecX.at[j][0], i) ; 
				count++; 
			}
			D.at[i][i] = Sx.at[2 * i][0];

			// 추가
			S.at[0][i] = Sx.at[i][0];
			S.at[i][0] = Sx.at[i][0];

		}

		S = addMat(S, D);

		
		for (int i = 1; i < n; i++) {

			S.at[n][i] = Sx.at[n + i][0];
			S.at[i][n] = Sx.at[n + i][0];

			for (int j = i + 1; j < n ; j++) { 
				
				S.at[i][j] = Sx.at[i + j][0]; 
				S.at[j][i] = Sx.at[i + j][0]; 
				count++;
			}
		}


		//printMat(S, " S in polyFit_mat_count ") ;

		
		invMat(S,S_inv);
		
		multMat_void(S_inv, b, vecZ);

	}

	printf("count : %d \n", count);

	return vecZ;
}






Matrix expFit(Matrix _vecX, Matrix _vecY) { 

	int n = 1;

	Matrix vecZ = createMat(n+1, 1);
	

	double a1 = 0.0;
	double a0 = 0.0;

	int mx = _vecX.rows;
	int my = _vecY.rows;


	if (mx != my) {
		printf("ERROR: The number of elements in x must be the same as in y.");
	}
	else 
	{
		vecZ = linearFit_mat(_vecX, _vecY) ; 
	}

	vecZ.at[0][0] = exp(vecZ.at[0][0]) ; 

	return vecZ ;
}


Matrix nonlinearSys(Matrix Funcs(Matrix _Z), Matrix Jacob(Matrix _Z), Matrix _Z0, double tol) {

	// Initialization
	int n = _Z0.rows;
	double loss = 10;

	Matrix J = zeros(n, n); 
	Matrix F = zeros(n, 1);  
	Matrix H = zeros(n, 1); 
	Matrix Z = zeros(n, 1); 

	Matrix d = zeros(n, 1); 
	Matrix U = eye(n, n); 

	// Initial condition
	copyMat(Z, _Z0);

	for (int k = 0; k < 20; k++) {


		// update F, J 
		F = Funcs(Z);
		J = Jacob(Z);

		printMat(J,"[ J ]");

		// F =-1*F;
		F = smultMat(F, -1.0);

		// Solving JH = -F using Gauss elimination
		gaussElim_Basic(J, F, U, d);
		backsub(U, d, H);

		// update Z
		//	Z = Z + H;
		Z = addMat(Z, H);

		// update F After Z update
		F = Funcs(Z);

		// calculate loss



		//loss = sqrt(loss);

		loss = norm_2(F);


		//printf("iter =%d \t x=%0.3f \t y=%0.3f \t loss=%0.3f\n", k, Z.at[0][0], Z.at[1][0], loss);

		if (loss < tol) {
			printf("Early Termination\n\r");
			return Z;
		}

	}
	return Z;
}
