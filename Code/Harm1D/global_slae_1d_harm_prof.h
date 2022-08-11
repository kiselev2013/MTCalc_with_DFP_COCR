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
 *  This file contains headers for the routines for SLAE and RHS assembly in 1D MT problem 
 *  Matrix is stored in skyline format 
 *                                                                                                                                          
 *  Written by Ph.D. Petr A. Domnikov                                                     
 *  Novosibirsk State Technical University,                                               
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                     
 *  p_domnikov@mail.ru                                                                    
 *  Version 1.2 November 21, 2020                                                         
*/                                                                                        


#pragma once
class Global_slae_1d_harm_prof
{
public:

	long n; // SLAE dimension
	long n_elem; // number of elements in the grid

	long* ig; // array of addresses of the beginning of lines
	double* ggl; // elements of the lower triangle
	double* ggu; // top triangle elements
	double* di; // diagonal
	double* pr; // right side vector

	long ig_n_1; // dimension ggl, ggu

	long* nvkat; // material numbers

	double* xyz; // grid

	long ay_hy;

	// equation coefficients
	double *mu;
	double *sigma;
	double omega;

	Global_slae_1d_harm_prof(long n, long n_elem, long *ig, long *nvkat,
		double *mu, double *sigma, double omega, double *xyz, long ay_hy);
	~Global_slae_1d_harm_prof();

	void Assembling_for_1d_harm_problem();
	void Set_boundary_conditions_for_Ay();
	void Set_boundary_conditions_for_Hy();
};
