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


#include "stdafx.h"
#include "portrait_harm_prof.h"
//-----------------------------------------------------------
Portrait_profil::Portrait_profil(long n_elem)
{
	this->n_elem = n_elem;
	this->mem = false;
}
//-----------------------------------------------------------
Portrait_profil::~Portrait_profil()
{
	if(this->mem == true)
	{
		delete [] this->ig;
	}
}
//-----------------------------------------------------------
void Portrait_profil::Gen_ig_for_1d_harm_problem()
{
	long i, k;

	this->n = (this->n_elem + 1)*2;

	this->ig = new long[n + 1];
	if(ig == 0)
	{
		char var[] = {"ig"};
		char func[] = {"Portrait_profil::Gen_ig_for_1d_harm_problem"};
		printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
		exit(1);
	}
	this->mem = true;

	ig[0] = 0;
	ig[1] = 0;
	ig[2] = 1;
	k = 3;

	for(i=0; i<this->n_elem; i++)
	{
		ig[k] = ig[k-1] + 2;
		ig[k+1] = ig[k] + 3;
		k += 2; 		
	}
}
//-----------------------------------------------------------