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
 *    Solving the SLAE in 1D MT by direct methods                                              
 *                                                                                             
 *  Written by Ph.D. Petr A. Domnikov                                                          
 *  Novosibirsk State Technical University,                                                    
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                          
 *  p_domnikov@mail.ru                                                                         
 *  Version 1.2 December 5, 2020                                                               
*/                                                                                             

#pragma once
class Solver_profil
{
public:
	long n;
	long *ig;
	double *di;
	double *ggl;
	double *ggu;
	double *pr;

	double *x;

    double zero; 

	double *h1, *h2; 
	bool mem_h; 

	Solver_profil();
	Solver_profil(long n, long *ig, double *di, double *ggl ,double *ggu, double *pr, double *x);
	~Solver_profil();


	int LU_profil();
	int Solve_SLAE_using_LU();
	int Solve_SLAE_using_LU(long n, long *ig, double *di, double *ggl ,double *ggu, double *pr, double *x);

	int LU_profil(long n, long *ig, double *di, double *ggl ,double *ggu, double *sl, double *su, double *d);
	int LLT_profil(long n, long *ig, double *di, double *gg, double *sl, double *d);

	double Scal(double *a, double *b, long n);

	void Solve_L_1(long n, long *ig, double *ggl, double *pr, double *y);

	int Solve_U(long n, long *ig, double *di, double *ggu, double *y, double *x, double *h);

	int Solve_L(long n, long *ig, double *di, double *ggl, double *y, double *x);

	int Solve_SLAE_using_LLT();

	int Solve_SLAE_using_LLT(long n, long *ig, double *di, double *ggl, double *pr, double *x);


};
