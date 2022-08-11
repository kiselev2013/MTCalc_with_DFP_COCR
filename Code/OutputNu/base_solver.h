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
 *  This file contains headers of subroutines for basic vector and matrix-vector operations   
 *                                                                                            
 *  Written by Ph.D. Petr A. Domnikov                                                         
 *  Novosibirsk State Technical University,                                                   
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                         
 *  p_domnikov@mail.ru                                                                        
 *  Version 1.3 April 7, 2021                                                                 
*/                                                                                            


#pragma once

class Base_solver
{
// The Base_solver class contains elementary routines for constructing iterative solvers (dot product, vector norm, etc.)          

protected:
	double *y_omp;


public:
	Base_solver();
	~Base_solver();

	int n_threads; // number of threads (processors for OpenMP)

	double Scal(double *x, double *y, long n);     // dot product                           
	double Norm_Euclid(double *x, long n);         // vector euclidean norm                 
	double Projection(double *vec, double *axis);  // projection of a vector onto an axis   
	double Relative_Error(double *analytic, double *numeric, long n); // relative error           
	double Spline(double x, long n, double *xyz, double *values);     // linear interpolation     

	// Matrix-vector multiplication  (in dense format)    
	void Mult_Plot(double *a,double *pr,double *rez,long n); 

	// Givens rotation   
	int  Givens1(double& x, double& y, double& c, double& s);
	void Givens2(double& x, double& y, double c, double s);
	int  Givens(double *a, double *f, long n);

	// Solution of a system with a lower triangular matrix in a dense format    
	int Undirect(double *a, double *b, double *x, long n);

	// Solution of a SLAE with a square matrix whose lower triangle contains only one non - zero subdiagonal.	 
	// This is necessary if the Arnoldi orthogonalization fails                                                      	
	int Solve_square_subdiag(double *a, double *b, double *x, long n);

	// writing to the file the relative discrepancy with which they came out, eps, the number of iterations and the time of solving the SLAE       
	int Write_kit(char *fname, double residual, double eps, long iter, __time64_t time);
	int Write_kit(char *fname, double residual, double eps, long iter, __time64_t time, long n);
	int Write_kit(char *fname, double residual, double eps, long iter, __time64_t time, long n, long size_jg);
	int Write_kit(char *fname, double residual, double eps, long iter, __time64_t time, double change_of_solution);


	// for complex-valued arithmetic (vectors are stored in the usual format, complex numbers (scalars) - in std::complex<double>)     
	                                                                                                                                   
	// dot product for complex-valued vectors                                                                                          
	std::complex<double> ScalCmplx(double *x, double *y, long nb); 

	// multiplication of a vector by a complex number          
	void MultCmplxNumVect(std::complex<double> a, double *x, double *y, long nb);

	// multiplication of components of one complex vector by components of another complex vector  
	void MultCmplxVectVect(long nb, double *a, double *b, double *c);

	// division of the components of one complex vector into components of another complex vector  
	void DivCmplxVectVect(long nb, double *a, double *b, double *c);

protected:
	long nMRSrestart; // number of iterations after which MRS is restarted
	// (This is to avoid rounding errors.)
public:
	void Set_nMRSrestart(long nMRSrestart) {this->nMRSrestart = nMRSrestart;}

};
//----------------------------------------------------
