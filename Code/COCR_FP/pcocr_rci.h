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
 *  This file contains the headers of the COCR solver routines.                                                   
 *  The reverse communication inteface is used. Matrix storage format and preconditioner is not specified here.   
 *                                                                                                             
 *  Written by Ph.D. Petr A. Domnikov                                                                          
 *  Novosibirsk State Technical University,                                                                    
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                          
 *  p_domnikov@mail.ru                                                                                         
 *  Version 1.2 April 7, 2021                                                                                  
*/                                                                                                             

#pragma once
#include "rci.h"
//------------------------------------------------------------------------
class PCOCR_RCI : public RCI
{
private:
	int nb;
	std::complex<double> alpha, betta;
	std::complex<double> zs, aw, zs_1;
	double *r; // initial system residual
	double *p; // search direction
	double *s; // residual of the preconditioned system
	double *z;  // z=A*s
	double *a;  // a=A*p
	double *w;  // w=M^{-1}*a

public:
	PCOCR_RCI(int n, int maxiter, double eps, double *x, double *pr, double **in, double **out);
	~PCOCR_RCI();

	int Run();
};
