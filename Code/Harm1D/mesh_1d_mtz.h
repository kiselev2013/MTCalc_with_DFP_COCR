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
 *  1D mesh generator                                                                     
 *                                                                                        
 *  Written by Ph.D. Petr A. Domnikov                                                     
 *  Novosibirsk State Technical University,                                               
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                     
 *  p_domnikov@mail.ru                                                                    
 *  Version 1.2 October 23, 2020                                                          
*/                                                                                        

#pragma once
const double X0_1D = -10000000.0; // lower boundary of the computational domain
const double X1_1D = 0.0;
const double COEF_RAZR = 1.02; // stretch factor
const double STEP0 = 1e-2;     // initial step size
//------------------------------------------------
class Mesh_1d_mtz
{
public:
	long n_points; 
	long n_elem; 
	long n_materials;
	long n_layers;
	double *coords;
	long *nvkat;
	double *layers;
	long *nodes_in_layer;

	double *mu;
	double *sigma;
	double omega;

	Mesh_1d_mtz();
	~Mesh_1d_mtz();

	int Read_1d_data_for_3d_problem();
	void Gen_1d_mesh();
};
