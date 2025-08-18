/* myNP_tutorial.h */

/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Ryu Min Seo
Created          : 08-30-2024
Modified         : 09-17-2024
Language/ver     : C in MSVS2019

Description      : myNP_tutorial.h
/----------------------------------------------------------------*/

#ifndef		_MY_NP_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NP_H
#define		PI		3.14159265358979323846264338327950288419716939937510582
#define     DEG2RAD(X) (X * PI / 180.0)




#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// print the vector
extern void printVec(double* vec, int row);

// calculates power x^N
extern double power(double _x, int N);

// calculates factorial x!
extern double factorial(int _x);



// returns sin(x) in unit[rad]
extern double sinTaylor(double _x);

// returns sin(x) in unit[deg]
extern double sindTaylor(double _x);

// returns cos(x) in unit[rad]
extern double cosTaylor(double _x);

//Taylor series approximation for exp(x)
extern double expTaylor(double _x);



// Non-LInear

// Bisection Method to solve non-linear eqation.
double bisection(double func(double n1),  double _a0, double _b0, double _tol);

// Newton-Raphson Method to solve non-linear eqation.
double newtonRaphson(double func(double n1), double dfunc(double n2), double _x0, double _tol);

// secant Method to solve non-linear eqation.
double secant(double func(double n), double _x0, double _x1, double _tol);



//Differentiation 


// 2-Point central  difference,3-Point forward / backward  difference
void gradientFunc(double func(const double x), double _x[], double dydx[], int m);

// 2-Point central  difference, 3-Point forward / backward  difference
void gradient1D(double _x[], double _y[], double dydx[], int m);

// Second Differentiation
void acceleration(double x[], double y[], double dy2dx2[], int m);



//Integration


// Simpson13 Method to integration with passing function (N is random)
double integral(double func(const double x), double a, double b, int n);

// Rectangluar Method to integration
double IntegrateRect(double x[], double y[], int m);

// Trapezoid Method to integration
double trapz(double x[], double y[], int m);

// Simpson13 Method to integration
double simpson13(double x[], double y[], int m);




//Interpolation

// Interpolation lagrange
double Linespline_lagrange(double x[], double y[], double xq , int m);

// Interpolation Newtonian
double Linespline_Newtonian(double x[], double y[], double xq, int m);

// Interpolation Standard 
double Linespline_Standard(double x[], double y[], double xq, int m);

//Linespline with Nst lagrange Polynomial method
double lagrange_spline(double x[], double y[], int m, int n, double xq);







// 기타 등등. 

//Q2 
double diff_a_p(double func(double n), double x, double h);

// Q_2 (2차미분)
double double_diff_a_p(double func(double n), double x, double h);



void odeEU(double y[], double odeFunc(const double t, const double y), const double t0, const double tf, const double h, const double y_init);



void odeRK2(double y[], double odeFunc(const double t, const double y),	const double t0, const double tf, const double h, const double y_init);


void odeRK3(double y[], double odeFunc(const double t, const double y), const double t0, const double tf, const double h, const double y0); 


void sys2RK2(double y1[], double y2[], void odeFuncSys(double dYdt[], const double t, const double Y[]), const double t0, const double tf, const double h, const double y1_init, const double y2_init);


void sys2RK4(double y1[], double y2[], void odeFuncSys(double dYdt[], const double t, const double Y[]), const double t0, const double tf, const double h, const double y1_init, const double y2_init);


//double newtonian(double x[], double y[], int m, int n, double xq); 

//void partial_x(double func(double x, double y), double** fx, double** fy, double** result, int m, int n);

#endif