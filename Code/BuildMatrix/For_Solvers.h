/**                                                                                       
 * GENERAL REMARKS                                                                        
 *                                                                                        
 *  This code is freely available under the following conditions:                         
 *                                                                                        
 *  1) The code is to be used only for non-commercial purposes.                           
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer. 
 *                                                                                        
 *  			MTCalc_with_DFP_COCR                                              
 *  This file contains some basic routines for vector-matrix operations and error messages
 *                                                                                        
 *  Written by Ph.D. Petr A. Domnikov                                                     
 *  Novosibirsk State Technical University,                                               
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                     
 *  p_domnikov@mail.ru                                                                    
 *  Version 1.2 April 7, 2021                                                             
*/                                                                                        

#pragma once
//------------------------------------------------------------------------    
double Scal(double *a, double *b, long n);  	// Dot product 
double Norm_Euclid(double *a, long n); 		// Euclidian norm of a vector 
double Norm_Max(double *a, long n);		// Max norm of a vector
double Projection_On_Axis(double *v,double *o);	//Projection of the vector v onto the o axis  
void Mult_Plot(double *a, double *x, double *y, long n);	// Matrix-vector multiplication in dense format  
void Mult_MV(long *ig, long *jg, double *ggl, double *ggu, double *di, double *x, double *y, long n);  // Matrix-vector multiplication in sparse format 
long Max_Long(long a, long b);	// Maximum of 2 numbers 
long Min_Long(long a, long b);	// Minimum of 2 numbers     
void Sort2(long *a, long *b);	// Sort 2 numbers
double Interval(double *x, double *y);	// Distance between two 3D-points
double Interval_Parallel_Lines(double *a0, double *a1, double *b0, double *b1);	// Distance between two parallel lines in 3D      
double Spline(double x, long n, double *xyz, double *values); // Linear interpolation 
double Calc_dof(double *J, double *func, long n_local_edge);
// reports a memory allocation error
void Memory_allocation_error(const char *var, const char *func);
// reports a file open error and throws an exception
void Cannot_open_file(const char *fname, const char *func);
// reports an error opening the file, but the program continues
void Cannot_open_file_but_continue(const char *fname, const char *func);
// for selecting sin-, cos-components from a non-stationary problem at one point     
int Spline_sin_cos(double *x, double *t, int m, double w, double *e);
void Mult_Plot_AV(double *a, double *x, double *y, long n, long m);
