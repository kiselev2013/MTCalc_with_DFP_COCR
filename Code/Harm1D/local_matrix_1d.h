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
 *  Calculations of local matrices and local RHS for 1d problem                            
 *                                                                                         
 *  Written by Ph.D. Petr A. Domnikov                                                      
 *  Novosibirsk State Technical University,                                                
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                      
 *  p_domnikov@mail.ru                                                                     
 *  Version 1.2 October 13, 2020                                                           
*/                                                                                         

#pragma once
class Local_matrix_1d
{
public:

	double b[2][2]; // stiffness matrix
	double c[2][2]; // mass matrix
	double f_re[2], f_im[2]; // right side values in nodes
	double g_re[2], g_im[2]; // right side vectors for Re and Im equations
	double a[4][4]; // matrix for the harmonic problem
	double g[4]; // vector of the right side for the harmonic problem

	// equation coefficients
	double mu;
	double sigma;
	double omega;

	double h; // size of finite element

	long alpha; // direction of the electric current


	Local_matrix_1d(double h, double mu, double sigma, double omega,
		double *f_re, double *f_im, long alpha);
	~Local_matrix_1d();

	void Compute_local_matrix_and_vector_harm();

	void Compute_local_matrix_harm();
	void Compute_local_vector_harm();
	void Compute_local_matrix_b();
	void Compute_local_matrix_c();
	void Compute_local_vector(double *f, double *g);
};
