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
 *  This file contains the routine to solve the 1D MT problem in the layered medium by FEM    
 *                                                                                              
 *  Written by Ph.D. Petr A. Domnikov                                                           
 *  Novosibirsk State Technical University,                                                     
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                           
 *  p_domnikov@mail.ru                                                                          
 *  Version 1.3 November 23, 2020                                                               
*/                                                                                              

#pragma once
class Divgrad_1d
{
public:
	double *coords_1d; 
	double *sin_1d;
	double *cos_1d;
	long n_1d; // number of nodes in a one-dimensional grid
	double bak; // lower tank limit
	double step0; // initial step (near the Earth's surface) for a one-dimensional grid
	double coef_razr; // discharging factor
	long n_layers_1d; // number of layers
	double* layers_1d; // layer boundaries	double *sigma_1d;
	long alpha;

	int Solve_1d_Problem_for_3d_task();
	
	Divgrad_1d();
	~Divgrad_1d();
};
