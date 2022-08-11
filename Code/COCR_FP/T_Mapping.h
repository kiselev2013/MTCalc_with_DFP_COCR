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
 *  This file contains headers of the class for transition matrix used to store a nonconforming mesh
 *                                                                                                         
 *                                                                                                         
 *  Written by Ph.D. Petr A. Domnikov                                                                      
 *  Novosibirsk State Technical University,                                                                
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                      
 *  p_domnikov@mail.ru                                                                                     
 *  Version 1.5 March 7, 2021                                                                              
*/                                                                                                         


#pragma once
//--------------------------------------------------------------------------------------------------------------
class T_Mapping_Vec
{
public: 
	T_Mapping_Vec(int (*nver)[14], double (*xyz)[3], int kuzlov, int kpar);
	~T_Mapping_Vec();

	// transition matrix in sparse column format
	int* ig_t;
	int* jg_t;
	double* gg_t;

	int* ig_s, * jg_s; // SIGMA structure (auxiliary structure for T-matrix)
	double* s_val; // values ​​that will be required to form all non-zero components of the matrix T

	int n_c; // number of continuous functions (continuous)
	int n_dc; // number of discontinuous functions (discontinuous)
	intn; // total functions n = n_c + n_dc

	int kuzlov; // number of nodes in the grid
	int kpar; // number of boxes
	int(*nver)[14]; // cells listed by their nodes + terminal nodes + element type
	double(*xyz)[3]; // vertex coordinates (all together)

	int(*edges)[2]; // edges defined by 2 vertices
	int(*ed)[25]; // cells listed by their edges + terminal edges + element type 
};
//--------------------------------------------------------------------------------------------------