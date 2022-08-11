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
 *  The RSF_solver class contains the basic routines for
 *  matrix operations in sparse string format          
 *            
 *                                                                                        
 *  Written by Ph.D. Petr A. Domnikov                                                     
 *  Novosibirsk State Technical University,                                               
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                     
 *  p_domnikov@mail.ru                                                                    
 *  Version 1.3 March 13, 2021                                                          
*/                                                                                        
                                                                                          



#pragma once
#include "base_solver.h"

class RSF_solver : public Base_solver
{
public:
	RSF_solver();
	~RSF_solver();


	/////////////////////   incomplete factorization for Compressed Row Storage   ///////////////////////////////

	int LLT(long *ig, long *jg, double *gg, double *di, double *sg, double *d, long n);

	int LLTd1(long *ig, long *jg, double *gg, double *di, double *sg, double *d, long n);

	int LU_sq(long *ig, long *jg, double *ggl, double *ggu, double *di,
		double *sl, double *su, double *d, long n);

	// 1 - for the matrix L
	int LU(long *ig, long *jg, double *ggl, double *ggu, double *di,
		double *sl, double *su, double *d, long n);



	/////////////////////  matrix-vector multiplications  ////////////////////////////

	// multiplication of a symmetric sparse matrix by a vector
	void mult_symmetr(long *ig, long *jg, double *gg, double *di, double *x, double *y, long n);

	// multiplication of a nonsymmetric sparse matrix by a vector
	void mult_mv(long *ig, long *jg, double *ggl, double *ggu, double *di, double *x, double *y, long n);
	void mult_mv_omp(long *ig, long *jg, double *ggl, double *ggu, double *di, double *x, double *y, long n, double *y_omp);

	// multiplication by an upper triangular matrix (not 1 on the diagonal)
	void mult_u_d(long *ig, long *jg, double *ggu, double *di, double *x, double *y, long n);

	// multiplication by a lower triangular matrix (not 1 on the diagonal)
	void mult_l_d(long *ig, long *jg, double *ggl, double *di, double *x, double *y, long n);


	/////////////////// SLAE solution with triangular matrices ///////////////////////

	// SLAE solution with a lower triangular matrix (not 1 on the diagonal)
	void solve_l_d(long *ig, long *jg, double *ggl, double *di, double *f, double *x, long n);

	// SLAE solution with an upper triangular matrix (not 1 on the diagonal)
	void solve_u_d(long *ig, long *jg, double *ggu, double *di, double *f, double *x, double *s, long n);


	//////////////////   Arnoldi orthogonalization   ///////////////////////////

	int Arnoldi(long n, long *ig, long *jg, double *di, double *ggl, double *ggu,
		double *v, double *h, double *w, long m, long *p);

	int Arnoldi_diag(long n, long *ig, long *jg, double *di, double *ggl, double *ggu,
		double *v, double *h, double *w, long m, long *p,
		double *d, double *help);

	int Arnoldi_lu_sq(long n, long *ig, long *jg, double *di, double *ggl, double *ggu,
		double *v, double *h, double *w, long m, long *p,
		double *d, double *help, double *sl, double *su);

	int Arnoldi_3diag(long n, long *ig, long *jg, double *di, double *ggl, double *ggu,
		double *v, double *h, double *w, long m, long *p,
		double *d, double *help, double *sl, double *su, long *ig_d, long *jg_d);





};