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
 *   Nonzero pattern for matrix in skyline format (1D MT)                                    
 *                                                                                           
 *  Written by Ph.D. Petr A. Domnikov                                                        
 *  Novosibirsk State Technical University,                                                  
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                        
 *  p_domnikov@mail.ru                                                                       
 *  Version 1.2 October 13, 2020                                                             
*/                                                                                           

#pragma once
class Portrait_profil
{
public:
	long n;       
	long n_elem;  

	long *ig; 
	bool mem;

	Portrait_profil(long n_elem);
	~Portrait_profil();

	void Gen_ig_for_1d_harm_problem();
};
