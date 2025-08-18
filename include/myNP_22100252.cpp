/* myNP_tutorial.cpp */
/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Ryu Min Seo 22100252
Created          : 08-30-2024
Modified         : 09-17-2024
Language/ver     : C in MSVS2019

Description      : myNP_tutorial.cpp
/----------------------------------------------------------------*/

#include "myNP_22100252.h"

// print the vector
extern void printVec(double* vec, int size)
{
	for (int i = 0; i < size; i++)
		printf("Vector[%d] = %.1f \n", i, vec[i]);
	printf("\n");
}



// factorial function
extern double factorial(int N)
{
	double y = 1.0;
	for (int k = 2; k <= N; k++)
		y = y * (double)k;

	return y;
}

// power fuction
extern double power(double _x, int N)
{
	double y = 1;

	for (int i = 0; i < N; i++)
	{
		y = y * _x;
	}


	return y;
	// [TODO] add your algorithm here
}



//  Taylor series approximation for sin(x) (input unit: [rad])
extern double sinTaylor(double _x)
{
	int N_max = 20;
	double S_N = 0;

	for (int k = 0; k < N_max; k++) {

		double fac = factorial(2 * k + 1);
		S_N = S_N + power(-1, k) * power(_x, 2 * k + 1) / fac;
	}

	return S_N;
}

// Taylor series approximation for sin(x) (input unit: [deg])
extern double sindTaylor(double _x)
{
	_x = _x * PI / 180;

	int N_max = 20;
	double S_N = 0;

	for (int k = 0; k < N_max; k++) {

		S_N = S_N + power(-1, k) * power(_x, 2 * k + 1) / factorial(2 * k + 1);
	}

	return S_N;

}

//  Taylor series approximation for cos(x) (input unit: [rad])
extern double cosTaylor(double _x)
{
	int N_max = 20; 
	double S_N = 0; 

	for (int k = 0; k < N_max; k++) {

		S_N = S_N + power(-1, k) * power(_x, 2 * k) / factorial(2 * k);
	}

	return S_N;
}

//Taylor series approximation for exp(x)
extern double expTaylor(double _x) {

	int N_max = 20 ;
	double S_N = 0 ;

	for (int k = 0; k < N_max; k++) {

		S_N = S_N + power(_x, k) / factorial(k);

	}

	return S_N;
}





// Non-linear


//  Bisection Method with Passing a Function 
/* Bisection Method
	func     : function parameter #1
	_a      : initial value #1
	_b      : initial value #2
	_tol   : tolerance

	condition : 1. solution within [_a,_b]
				2. f(_a)*f(_b) < 0 
*/
double bisection(double func(double n1), double _a0, double _b0, double _tol) 
{
	// Define Bisection function, assuming (func(a) * func(b) <0 )
	// Initialization
	int k = 0 ; 
	int Nmax = 100 ; 
	double a = _a0 ; 
	double b = _b0 ; 
	double xn = 0.00  ; 
	double ep = 1000.00 ; 


	// Bisection 
	while (k<Nmax && ep>_tol) {

		// Update xn as midpoint
		xn = (a + b) / 2;

		if (func(a)*func(b)>0) {
			return 0.00 ;
		}

		// Update range a, b
		if (func(xn) * func(a) < 0) 
		{
			b = xn ;
		}
		else if (func(xn) * func(a) > 0) 
		{
			a = xn ;
		}

		// Check tolerance
		ep = fabs(func(xn)); 

		k++;

		printf("k:%d \t", k);
		printf("Xn(k): %f \t", xn);
		printf("Tol: %.8f\n", ep);
	}

	return xn ;
}


// Newton-Raphson Method with Passing a Function 
/* Newton-Raphson
	func     : function parameter #1
	dfunc    : function parameter #2
	_x0      : initial value 
	_tol     : tolerance

	condition : 
	1. initial point(_x0) have to be  nearby true value
	2. dfubc(xi) != 0

*/
double newtonRaphson(double func (double n1), double dfunc(double n2), double _x0, double _tol)
{
	double xn = _x0 ;
	double ep = 1000;
	int Nmax = 1000;
	int k = 0;
	double h = 0;

	while (k<Nmax && ep>_tol) {
		if (dfunc(xn) == 0)
		{
			printf("[ERROR] dF == 0 !!\n");
			break;
		}
		else
		{
			// get h
			h = -1 * func(xn) / dfunc(xn) ; 

			// update x(k+1)=x(k)+h(k)
			xn = xn + h ; 

			// check tolerance 
			ep = fabs(func(xn)) ;

			k++;

			printf("k:%d \t", k);
			printf("X(k): %f \t", xn);
			printf("Tol: %.10f\n", ep);

		}
	}

	// Case of failure to find solution until k has become Nmax. 
	if (k == Nmax) {
		printf("A solution could not be found until k has become Nmax.\n");
	}

	return xn ;
}


// Secant Method with Passing a Function 
/* Secant
	func     : function parameter 
	_x0      : initial value #1
	_x1      : initial value #2
	_tol     : tolerance
*/
double secant(double func(double n), double _x0, double _x1, double _tol)
{
	double xn_1 = _x0 ;
	double xn = _x1 ;
	double xn1 = 0;
	double ep = 1000;
	int Nmax = 1000;
	int k = 0;

	while (k<Nmax && ep>_tol) {
		if (func(xn) - func(xn_1) == 0) // Error check
		{
			printf("[ERROR] dF == 0 !!\n");
			break;
		}
		else
		{
			xn1 = xn - func(xn) * ((xn - xn_1) / (func(xn) - func(xn_1)));	// Update xn1
			ep = fabs(func(xn1));   // Check tolerance 

			xn_1 = xn; // Update xn_1
			xn = xn1;  // Update xn, 
			k++;       // Update k

			printf("k:%d \t", k);
			printf("X(k): %f \t", xn);
			printf("Tol: %.10f\n", ep);
		}
	}
	if (k == Nmax) // Case of failure to find solution until k has become Nmax.
	{
		printf("A solution could not be found until k has become Nmax.\n");
	}


	return xn ;
}





// differentiation

// 2-Point central  difference, 3-Point forward / backward  difference  (1차 미분) with passing function.
/* o(h)2 method
	func     : passing function
	x[]      : data set of x
	dxdy[]   : data set of diffence value
	m        : # of datas
*/
void gradientFunc(double func(const double x), double x[], double dydx[], int m) {

	dydx[0] = (-3 * func(x[0]) + 4 * func(x[1]) - func(x[2])) / (x[2] - x[0]);

	for (int i = 1; i < m - 1; i++)
	{
		dydx[i] = (func(x[i + 1]) - func(x[i - 1])) / (x[i + 1] - x[i - 1]);
	}

	dydx[m - 1] = (3 * func(x[m - 1]) - 4 * func(x[m - 2]) + func(x[m - 3])) / (x[m - 1] - x[m - 3]);

}





// 2-Point central  difference, 3-Point forward / backward  difference (1차 미분)
/* o(h)2 method
	x[]      : data set of x
	y[]      : data set of y
	dxdy[]   : data set of diffence value
	m        : # of datas
*/
void gradient1D(double x[], double y[], double dydx[], int m) {

	// 3-Point forward difference
	dydx[0] = (-3 * y[0] + 4 * y[1] - y[2]) / (x[2] - x[0]);

	// 2-Point central  difference
	for (int i = 1; i < m - 1; i++)
	{
		dydx[i] = (y[i+1] - y[i - 1]) / (x[i + 1] - x[i-1]);
	}

	// 3-Point backward  difference
	dydx[m-1] = (3 * y[m-1] - 4 * y[m-2] + y[m-3]) / (x[m-1] - x[m-3]);
}

// 3-Point central  difference, 4-Point forward / backward  difference (2차 미분)
/* Second Differentiation 
	x[]      : data set of x
	y[]      : data set of y
	dx2dy2[]   : data set of diffence value
	m        : # of datas
*/
void acceleration(double x[], double y[], double dy2dx2[], int m) {

	// 4-Point forward difference
	dy2dx2[0] = (2 * y[0] - 5 * y[1] + 4 * y[2] - y[3]) / pow((x[1] - x[0]),2);

	// 3-Point central  difference
	for (int i = 1; i < m - 1; i++)
	{
		dy2dx2[i] = (y[i+1] - 2* y[i] + y[i-1]) / pow((x[i+1] - x[i]),2);
	}

	// 4-Point backward  difference
	dy2dx2[m - 1] = (2 * y[m-1] - 5 * y[m-2] + 4 * y[m-3] - y[m-4]) / pow((x[m-1] - x[m-2]), 2);
}





// Integration

// Simpson13 + Simpson38  Method to integration with passing function
/* Simpson13 + Simpson38Method
	func     : passing function
	a        : starting point
	a        : end point
	n        : # of datas
*/
double integral(double func(const double x), double a, double b, int n) {

	int N = n - 1;  // interval 개수 설정
	double I = 0;  // 적분 받을 변수
	double h = (b - a) / N; // 간격 h 설정

	// 1/3 simpson method

	// interval 개수가 짝수일 경우 (1/3 simpson method)
	if (N % 2 == 0) {

		double f2i = 0; // 2번 반복되는f(x2i) 합해주는 변수

		for (int i = 1; i < N / 2; i++) {
			f2i = f2i + func(a + 2 * i * h);
		}

		double f2i_1 = 0; // 4로 곱해지는 f(x2i) 합해주는 변수

		for (int i = 1; i < (N / 2) + 1; i++) {
			f2i_1 = f2i_1 + func((2 * i * h) - 1 * h);
		}

		I = h / 3 * (func(a) + func(b) + 2 * f2i + 4 * f2i_1);

	}

	// interval 개수가 홀수일 경우 (1/3 simpson method + 3/8 simpson method)
	else {

		// 1/3 simpson method

		double f2i = 0; // 2번 반복되는f(x2i) 합해주는 변수

		for (int i = 1; i < ((N - 1) / 2) - 1; i++) {
			f2i = f2i + func(a + 2 * i * h);
		}


		double f2i_1 = 0; // 4로 곱해지는 f(x2i) 합해주는 변수

		for (int i = 1; i < ((N - 1) / 2); i++) {
			f2i_1 = f2i_1 + func(a + (2 * i * h) - 1 * h);
		}

		I = h / 3 * (func(a) + func(b - 3 * h) + 2 * f2i + 4 * f2i_1);

		// 3/8 simpson method at last 4 points
		I = I + (3 * h * (func(b - 3 * h) + 3 * func(b - 2 * h) + 3 * func(b - h) + func(b))) / 8;

	}

	return I;

}



// Rectangluar Method to integration
/* Rectangluar Method
	x[]      : data set of x
	y[]      : data set of y
	m        : # of datas
*/
double IntegrateRect(double x[], double y[], int m) {

	int N = m - 1;
	double I = 0;

	for (int i = 0; i < N; i++)
		I += y[i] * (x[i + 1] - x[i]) ;
	return I;
}

// Trapezoid Method to integration
/* Trapezoid Method
	x[]      : data set of x
	y[]      : data set of y
	m        : # of datas
*/
double trapz(double x[], double y[], int m) { 

	int N = m - 1;
	double I = 0; 

	for (int i = 0; i < N; i++) {
		I = I + ((y[i] + y[i + 1]) / 2) * (x[i + 1] - x[i]);
	}

	return I;
}


// Simpson13 Method to integration
/* Simpson13 Method
	x[]      : data set of x
	y[]      : data set of y
	m        : # of datas
*/
double simpson13(double x[], double y[], int m) {

	int N = m - 1; // 인터벌 개수

	if (N % 2 != 0) {
		printf("Interval 개수(N)가 홀수입니다. 0으로 리턴합니다.  \n");
		return 0.0 ;
	}

	double I = 0;
	double h = x[1] - x[0];

	double f2i = 0;

	for (int i = 1; i < N / 2 ; i++) {
		f2i = f2i + y[2 * i];
	}

	double f2i_1 = 0;

	for (int i = 1; i < (N / 2) + 1; i++) {
		f2i_1 = f2i_1 + y[(2 * i) - 1];
	}

	I = h / 3 * (y[0] + y[N] + 2 * f2i + 4 * f2i_1);
	
	return I;

}








// Interpolation



//Linespline with lagrange Polynomial method
double Linespline_lagrange(double x[], double y[], double xq, int m) {

	int k = 0;
	double result = 0;
	
	for (int i = 0; i < m-1; i++) {
		
		if (x[i] <= xq && xq <= x[i + 1]) {
			k = i;
			
		}
	}

	printf("%d\n", k);
	result = y[k] * ((xq - x[k + 1]) / (x[k] - x[k + 1])) + y[k + 1] * ((xq - x[k]) / (x[k + 1] - x[k]));
	

	return result;

}



//Linespline with Newton polymiial method
double Linespline_Newtonian(double x[], double y[], double xq, int m) {

	int k = 0;

	int a = 0;
	for (int i = 0; i < m - 1; i++) {

		//a = (x[i] <= xq && xq <= x[i + 1]);
		if (x[i] <= xq && xq <= x[i + 1]) {
			k = i;
			printf("%d\n", k);
		}

	}

	double result = y[k] + (xq - x[k]) * ((y[k + 1] - y[k]) / (x[k + 1] - x[k]));

	return result;

}


//Linespline with standard polymiial method
double Linespline_Standard(double x[], double y[], double xq, int m) {

	int k = 0;

	for (int i = 0; i < m - 1; i++) {

		if (x[i] <= xq && xq < x[i + 1]) {
			k = i;
			printf("%d", k);
		}

	}

	double result = ((y[k + 1] - y[k]) / (x[k + 1] - x[k])) * xq - ((y[k + 1] - y[k]) / (x[k + 1] - x[k])) * x[k] + y[k] ; 

	return result;
}


//Linespline with Nst lagrange Polynomial method
double lagrange_spline(double x[], double y[], int m, int n, double xq) {

	int k = 0;

	for (int i = 0; i < m; i = i + n) {
		if (x[i] <= xq && xq <= x[i + n]) {
			k = i;
		}
	}

	// n = n차 방정식
	double yq = 0.0;
	double a = 1.0;

	printf("%d\n", k);

	for (int i = 0; i <= n; i++) {

		a = y[k + i];

		for (int j = 0; j <= n; j++) {
			if (j == i) continue;

			a = a * (xq - x[k + j]) / (x[k + i] - x[k + j]);
		}

		yq = yq + a;
	}

	return yq;
}



// Q_2 (1차 미분) 
double diff_a_p(double func(double n), double x, double h) {

	double y = (func(x + h) - func(x - h)) / (2 * h) ;

	return y;

}

// Q_2 (2차미분)
double double_diff_a_p(double func(double n), double x, double h) {
	
	double y = (func(x + h) + 2 * func(x) + func(x - h)) / pow(h,2) ;

	return y ; 

}




void odeEU(double y[], double odeFunc(const double t, const double y), const double t0, const double tf, const double h, const double y_init)
{
	// 
	// [BRIEF DESCRIPTION OF THE FUNCTION  GOES HERE]
	// 
		// Variable Initialization
	int N = (tf - t0) / h + 1;

	//Initial Condition
	double ti = t0;
	y[0] = y_init;

	// Euler Explicit ODE Method
	for (int i = 0; i < N - 1; i++) {
		
		y[i + 1] = y[i] + odeFunc(ti, y[i]) * h ;
		ti += h ;
	}
}


void odeRK2(double y[], double odeFunc(const double t, const double y), const double t0, const double tf, const double h, const double y_init)
{


	int N = (tf - t0) / h + 1;
	double y2 = 0;
	double t2 = 0;

	// Initial Condition
	double ti = t0;
	y[0] = y_init;

	// RK Design Parameters
	double C1 = 0.5;
	double C2 = 1 - C1;
	double alpha = 1;
	double beta = alpha;
	double K1 = 0 ;
	double K2 = 0 ;

	// RK2 ODE Solver
	for (int i = 0; i < N - 1; i++) {

		// [First-point Gradient]
		K1 = odeFunc(ti, y[i]);

		// [Second-point Gradient]
		t2 = ti + alpha * h ;
		y2 = y[i] + beta* K1 * h ;

		K2 = odeFunc(t2, y2); 

		y[i+1] = y[i] + (C1 * K1 + C2 * K2) * h ;

		ti += h ;

	}
}


void odeRK3(double y[], double odeFunc(const double t, const double y), const double t0, const double tf, const double h, const double y0)
{
	// 
	// [BRIEF DESCRIPTION OF THE FUNCTION  GOES HERE]
	// 

	int N = (tf - t0) / h + 1;
	double y2 = 0;
	double t2 = 0;

	double y3 = 0;
	double t3 = 0;

	// Initial Condition
	double ti = t0;
	y[0] = y0;

	// RK Design Parameters
	double C1 = 1.0 ;
	double C2 = 4.0 ;
	double C3 = 1.0 ;

	double alpha_1 = 0.5 ;
	double alpha_2 =  1.0  ;

	double beta_1 = 0.5 ;
	double beta_2 = - 1.0 ;
	double beta_3 = 2.0  ;

	double K1 = 0.0 ; 
	double K2 = 0.0 ; 
	double K3 = 0.0 ; 

	// RK2 ODE Solver
	for (int i = 0; i < N - 1; i++) {

		// [First-point Gradient]
		K1 = odeFunc(ti, y[i]); 

		t2 = ti + alpha_1 * h ; 
		y2 = y[i] + beta_1 * K1 * h ;

		K2 = odeFunc(t2, y2);

		t3 = ti + alpha_2 * h;
		y3 = y[i] + beta_2 * K1 * h + beta_3 * K2 * h ;

		K3 = odeFunc(t3, y3);

		y[i + 1] = y[i] + (C1 * K1 + C2 * K2 + C3 * K3 ) * h / 6.0 ;

		ti += h ;

	}

}









// ODE RK2 for 2nd order ODE 
void sys2RK2(double y1[], double y2[], void odeFuncSys(double dYdt[], const double t, const double Y[]), const double t0, const double tf, const double h, const double y1_init, const double y2_init)
{


		// Variable Initialization
	int N = (tf - t0) / h + 1 ;
	double ti = t0 ;

	// Initial Condition	
	y1[0] = y1_init;  // Y
	y2[0] = y2_init;  // v


	double dYdt[2] = { 0 };
	double Yin[2] = { y1[0], y2[0] } ; 

	double Ky1 = 0.0;
	double Ky2 = 0.0;

	double Kv1 = 0.0;
	double Kv2 = 0.0;

	double t2 = 0;

	for (int i = 0; i < N - 1; i++) {

		Yin[0] = y1[i];
		Yin[1] = y2[i];

		odeFuncSys(dYdt, ti, Yin);

		Ky1 = dYdt[0];
		Kv1 = dYdt[1];


		Yin[0] = y1[i] + Ky1 * h;
		Yin[1] = y2[i] + Kv1 * h;

		t2 = ti + h;

		odeFuncSys(dYdt, t2, Yin);

		Ky2 = dYdt[0];
		Kv2 = dYdt[1];


		y1[i + 1] = y1[i] + 0.5 * (Ky1 + Ky2) * h;
		y2[i + 1] = y2[i] + 0.5 * (Kv1 + Kv2) * h;

		ti = ti + h;

	}



	// RK2 ODE for 2nd ODE
	// [YOUR CODE GOES HERE]	
	// [YOUR CODE GOES HERE]	


}



void sys2RK4(double y1[], double y2[], void odeFuncSys(double dYdt[], const double t, const double Y[]), const double t0, const double tf, const double h, const double y1_init, const double y2_init) {


	int N = (tf - t0) / h + 1 ;
	double ti = t0 ;

	y1[0] = y1_init;  // Y
	y2[0] = y2_init;  // v


	double dYdt[2] = { 0 } ; 

	double Yin[2] = { 0 } ; 

	double Ky1 = 0.0 ;
	double Ky2 = 0.0 ;
	double Ky3 = 0.0 ;
	double Ky4 = 0.0 ;

	double Kv1 = 0.0 ;  
	double Kv2 = 0.0 ;
	double Kv3 = 0.0 ;
	double Kv4 = 0.0 ;
 
	double C1 = 1.0 / 6.0 ;
	double C2 = 2.0 / 6.0 ;
	double C3 = 2.0 / 6.0 ;
	double C4 = 1.0 / 6.0 ;

	double alpha_1 = 0.5 ;
	double alpha_2 = 0.5 ;
	double alpha_3 = 1.0 ;

	double beta_1 = 0.5 ;
	double beta_2 = 0.0 ;
	double beta_3 = 0.5 ;
	double beta_4 = 0.0 ;
	double beta_5 = 0.0 ;
	double beta_6 = 1.0 ;

	double t2 = 0;
	double t3 = 0;
	double t4 = 0;

	for (int i = 0; i < N - 1; i++) {


		//  t1 
		Yin[0] = y1[i] ;
		Yin[1] = y2[i] ;

		odeFuncSys(dYdt, ti, Yin) ;

		Ky1 = dYdt[0] ;
		Kv1 = dYdt[1] ;
		

		//  t2

		t2 = ti + alpha_1 * h ;

		Yin[0] = y1[i] + beta_1 * Ky1 * h ;
		Yin[1] = y2[i] + beta_1 * Kv1 * h ;

		odeFuncSys(dYdt, t2, Yin) ;

		Ky2 = dYdt[0] ;
		Kv2 = dYdt[1] ;

		//  t3

		t3 = ti + alpha_2 * h;  

		Yin[0] = y1[i] + beta_2 * Ky1 * h + beta_3 * Ky2 *  h;
		Yin[1] = y2[i] + beta_2 * Kv1 * h + beta_3 * Kv2 * h ;
		 
		odeFuncSys(dYdt, t3, Yin);  

		Ky3 = dYdt[0]; 
		Kv3 = dYdt[1];

		//  t4

		t4 = ti + alpha_3 * h;

		Yin[0] = y1[i] + beta_4 * Ky1 * h  + beta_5 * Ky2 * h  + beta_6 * Ky3 * h ;
		Yin[1] = y2[i] + beta_4 * Kv1 * h  + beta_5 * Kv2 * h  + beta_6 * Kv3 * h ;

		odeFuncSys(dYdt, t4, Yin);

		Ky4 = dYdt[0];
		Kv4 = dYdt[1];



		y1[i + 1] = y1[i] + (C1 * Ky1 + C2 * Ky2 + C3 * Ky3 + C4 * Ky4) * h;
		y2[i + 1] = y2[i] + (C1 * Kv1 + C2 * Kv2 + C3 * Kv3 + C4 * Kv4) * h;

		ti = ti + h;


	}


}



/* 
double newtonian(double x[], double y[], int m, int n, double xq) {
	int k = 0;

	for (int i = 0; i < m; i = i + n) {
		if (x[i] <= xq && xq <= x[i + n]) {
			k = i;
		}
	}

	// n = n차 방정식
	double yq = 0.0 ;
	double a = 1.0 ;

	yq = y[k];

	for (int i = 1; i <= n; i++) {

		for (int j = 0; j < i; j++) {
			a = a * (xq - x[k + j]);
		}


	}
		return yq;
}
*/



/*
void partial_x(double func(double x, double y) , double** fx, double** fy, double** result, int m, int n) {

	for (int i = 0 ; i < n ; i++) {
	
		for (int j = 0 ; i < m; j++) {
			if (j!=m-1) {
				result[i][0] = (func(fx[i][j+1], fy[i][j])- func(fx[i][j], fy[i][j])) / (fx[i][j+1]- fx[i][j]);
			}
			else {
				result[i][9] = (func(fx[i][9], fy[i][9]) - func(fx[i][8], fy[i][9])) / (fx[i][9] - fx[i][8]);
			}			
		}
	}
}
*/