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
 *  This file contains some basic routines for Solver: locally optimal scheme;
 *  matrix storage format: sparse row-column format                                         
 *                                                                                                 
 *  Written by Ph.D. Petr A. Domnikov                                                              
 *  Novosibirsk State Technical University,                                                        
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                              
 *  p_domnikov@mail.ru                                                                             
 *  Version 1.3 January 10, 2021                                                                  
*/                                                                                                 


#pragma once
#include "rsf_solver.h"

class LOS_rsf_solver : public RSF_solver
{
public:

	int ShowProgress;

	LOS_rsf_solver();
	~LOS_rsf_solver();

	long LOS_LU_sq(long n, long *ig, long *jg, double *di, double *ggl, double *ggu, double *pr,
		double *x, double eps, long maxiter);
	long LOS_LU_sq(long n, long *ig, long *jg, double *di,
		double *ggl, double *ggu, double *f, double *x,
		double eps, long maxiter, double *d, double *sl, double *su, 
		double *r, double *z, double *p, double *qr, double *laqr, double *h);

	long LOS_diag(long n, long *ig, long *jg, double *di, double *ggl, double *ggu, double *pr,
		double *x, double eps, long maxiter);
	long LOS_diag(long n, long *ig, long *jg, double *di,
		double *ggl, double *ggu, double *f, double *x,
		double eps, long maxiter, double *d, 
		double *r, double *z, double *p, double *qr, double *laqr, double *h, double *temp_x);

};