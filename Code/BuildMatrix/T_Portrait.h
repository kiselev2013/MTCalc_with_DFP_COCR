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
 *  This file contains some basic routines for construction of a matrix portrait of a finite element SLAE
 * (vector FEM for solving 3D problems)                                                                  
 *                                                                                                       
 *  Written by Ph.D. Petr A. Domnikov                                                                    
 *  Novosibirsk State Technical University,                                                              
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                    
 *  p_domnikov@mail.ru                                                                                   
 *  Version 1.3 April 6, 2021                                                                            
*/                                                                                                       


#pragma once

class T_Portrait
{
public:
	long (*ed)[25];
	long n_elem,n,size_jg;
	long *ig,*jg,*idi,*ijg; 
	T_Portrait(long *ed,long n,long n_elem);
	~T_Portrait();
	void Gen_Portrait();
	void Gen_idi_ijg(long *nvkat, long (*nver)[14]);
	void Set_type_of_block(long *target_array, long adr, long type);
};

const int FILTER_MASS_MATRIX_VEC[12][12] = { 
// 1  2  3  4  5  6  7  8  9  10 11 12
   2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
   2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
   2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
   2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 
   1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 
   1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 
   1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,  
   1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,  
   1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,  
   1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2
};
