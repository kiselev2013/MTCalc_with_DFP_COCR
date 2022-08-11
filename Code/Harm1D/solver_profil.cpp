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


#include "stdafx.h"
#include "solver_profil.h"
//-----------------------------------------------------------
Solver_profil::Solver_profil(long n, long *ig, double *di, double *ggl,
							 double *ggu, double *pr, double *x)
{
	this->n = n;
	this->ig = ig;
	this->di = di;
	this->ggl = ggl;
	this->ggu = ggu;
	this->pr = pr;
	this->x = x;

	this->zero = 1e-100;
	this->mem_h = false;
}
//-----------------------------------------------------------
Solver_profil::Solver_profil()
{
	this->zero = 1e-100;
	this->mem_h = false;
}
//-----------------------------------------------------------
Solver_profil::~Solver_profil()
{
	if(this->mem_h == true)
	{
		delete [] this->h1;
		delete [] this->h2;
	}
}
//-----------------------------------------------------------
int Solver_profil::LU_profil()
{
	return this->LU_profil(n, ig, di, ggl, ggu, ggl, ggu, di);
}
//-----------------------------------------------------------------
//---- SLAE solution (i.e. LU decomposition and then SLAE solution with ---
//---- triangular matrices L and U) ---
//-----------------------------------------------------------------
int Solver_profil::Solve_SLAE_using_LU()
{
	return this->Solve_SLAE_using_LU(n, ig, di, ggl, ggu, pr, x);
}
//-----------------------------------------------------------
int Solver_profil::Solve_SLAE_using_LU(long n, long *ig, double *di, double *ggl ,double *ggu, double *pr, double *x)
{
	if(this->LU_profil(n, ig, di, ggl, ggu, ggl, ggu, di)!=0)
		return 1;

	if(this->mem_h == false)
	{
		this->h1 = new double[n];
		if(h1 == 0)
		{
			char var[] = {"h1"};
			char func[] = {"Solver_profil::Solve_SLAE_using_LU"};
			printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
			exit(1);
		}

		this->h2 = new double[n];
		if(h2 == 0)
		{
			char var[] = {"h2"};
			char func[] = {"Solver_profil::Solve_SLAE_using_LU"};
			printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
			exit(1);
		}

		this->mem_h = true;
	}

	this->Solve_L_1(n, ig, ggl, pr, h1);

	if(this->Solve_U(n, ig, di, ggu, h1, x, h2)!=0)
		return 2;

	return 0;
}
//-----------------------------------------------------------
int Solver_profil::Solve_SLAE_using_LLT(long n, long *ig, double *di, double *ggl, double *pr, double *x)
{
	if(LLT_profil(n, ig, di, ggl, ggl, di)!=0)
		return 1;

	if(this->mem_h == false)
	{
		this->h1 = new double[n];
		if(h1 == 0)
		{
			char var[] = {"h1"};
			char func[] = {"Solver_profil::Solve_SLAE_using_LLT"};
			printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
			exit(1);
		}

		this->h2 = new double[n];
		if(h2 == 0)
		{
			char var[] = {"h2"};
			char func[] = {"Solver_profil::Solve_SLAE_using_LLT"};
			printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
			exit(1);
		}

		this->mem_h = true;
	}

	if(this->Solve_L(n, ig, di, ggl, pr, h1)!=0)
		return 2;

	if(this->Solve_U(n, ig, di, ggl, h1, x, h2)!=0)
		return 3;

	return 0;
}
//-----------------------------------------------------------
//------------ dot product  ----------------------
//-----------------------------------------------------------
double Solver_profil::Scal(double *x, double *y, long n)
{
	double sum = 0.0;
	long i;

	for(i=0;i<n;i++) sum += x[i]*y[i];

	return sum;
}
//-------------------------------------------------------------
//------   LU decomposition of a matrix in the skyline format
//-------------------------------------------------------------
int Solver_profil::LU_profil(long n, long *ig, double *di, double *ggl ,double *ggu, double *sl, double *su, double *d)
{
	long i, j, m;
	long razn; 
	long nonzero_i; 
	long n_elem_str_i, n_elem_str_j; 
	long str_beg_i, str_beg_j; 
	long str_end_i, str_end_j; 

	printf("Direct LU factorization... ");

	for(i=0; i<n; i++)
	{
		str_beg_i = ig[i];
		str_end_i = ig[i+1]-1;
		n_elem_str_i = str_end_i - str_beg_i + 1;

		for(m=str_beg_i; m<=str_end_i; m++)
		{
			j = i - n_elem_str_i + (m-str_beg_i);
			str_beg_j = ig[j];
			str_end_j = ig[j+1]-1;
			n_elem_str_j = str_end_j - str_beg_j + 1;
			nonzero_i = m - str_beg_i;

			razn = nonzero_i - n_elem_str_j;

			if(fabs(d[j]) < this->zero)
			{
				printf("Error: LU failed!!! Divizion by zero!!!\n");
				printf("i=%ld, j=%ld\n", i, j);
				return 1;
			}

			if(razn >= 0)
			{
				su[m] = ggu[m] - Scal(&sl[str_beg_j], &su[str_beg_i+razn], n_elem_str_j);	

				sl[m] = (ggl[m] - Scal(&sl[str_beg_i+razn], &su[str_beg_j], n_elem_str_j))/d[j];
			}
			else 
			{
				su[m] = ggu[m] - Scal(&sl[str_beg_j-razn], &su[str_beg_i], nonzero_i);	

				sl[m] = (ggl[m] - Scal(&sl[str_beg_i], &su[str_beg_j-razn], nonzero_i))/d[j];
			}
		}

		d[i] = (di[i] - Scal(&sl[str_beg_i], &su[str_beg_i], n_elem_str_i));
	}

	printf("done.\n");
	return 0;
}
//---------------------------------------------------------------
//-- SLAE solution with a lower triangular matrix (1 on the diagonal)
//---------------------------------------------------------------
void Solver_profil::Solve_L_1(long n, long *ig, double *ggl, double *pr, double *y)
{
	long i;
	long str_beg, n_elem_str;
	double s;

	for(i=0; i<n; i++)
	{
		str_beg = ig[i];
		n_elem_str = ig[i+1] - str_beg;

		s = Scal(&ggl[str_beg], &y[i-n_elem_str], n_elem_str);
        y[i] = (pr[i] - s);
	}
}
//-----------------------------------------------
//--SLAE solution with upper triangular matrix --
//-----------------------------------------------
int Solver_profil::Solve_U(long n, long *ig, double *di, double *ggu,
						   double *y, double *x, double *h)
{
	long i, j, k;

	for(i=0; i<n; i++) 
		h[i] = 0.0;

	for(i=n-1; i>=0; i--)
	{
		if(fabs(di[i]) < this->zero)
		{
			printf("Error: Solve_U failed!!! Divizion by zero!!!\n");
			return 1;
		}

		x[i] = (y[i] - h[i])/di[i];
		for(j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = i - ig[i+1] + j; 
			h[k] += ggu[j]*x[i];
		}
	}
	return 0;
}
//---------------------------------------------------------------
//---------------  Cholesky decomposition  -----------------------
//---------------------------------------------------------------
int Solver_profil::LLT_profil(long n, long *ig, double *di, double *gg, double *sl, double *d)
{
	long i, j, m;
	long str_beg_i, str_beg_j;       
	long str_end_i, str_end_j;       
	long n_elem_str_i, n_elem_str_j; 
	long nonzero_i; 
	long razn; 
	double s;

	printf("Direct LLT factorization... ");

	for(i=0; i<n; i++)
	{
		str_beg_i = ig[i];
		str_end_i = ig[i+1]-1;
		n_elem_str_i = str_end_i - str_beg_i + 1;

		j = i - n_elem_str_i; 
		for(m=str_beg_i; m<=str_end_i; m++)
		{	
			str_beg_j = ig[j];
			str_end_j = ig[j+1]-1;
			n_elem_str_j = str_end_j - str_beg_j + 1;
			nonzero_i = m - str_beg_i;

			if(fabs(d[j]) < zero)
			{
				printf("Error: LLT failed!!! Divizion by zero!!!\n");
				printf("i=%ld, j=%ld\n", i, j);
				return 1;
			}

			razn = nonzero_i - n_elem_str_j;

			if(razn >= 0)
			{
				sl[m] = (gg[m] - Scal(&sl[str_beg_i+razn], &sl[str_beg_j], n_elem_str_j))/d[j];
			}
			else 
			{
				sl[m] = (gg[m] - Scal(&sl[str_beg_i], &sl[str_beg_j-razn], nonzero_i))/d[j];
			}

			j++;
		}

		s = di[i] - Scal(&sl[str_beg_i], &sl[str_beg_i], n_elem_str_i);

		if(s <= 0.0)
		{
			printf("Error: LLT failed!!! Sqrt!!!\n");
			printf("i=%ld, j=%ld\n", i, j);
			return 1;
		}

		d[i] = sqrt(s);
	}
	printf("done.\n");
	return 0;
}
//---------------------------------------------------------------
//-- SLAE solution with a lower triangular matrix -------------------
//---------------------------------------------------------------
int Solver_profil::Solve_L(long n, long *ig, double *di, double *ggl, double *pr, double *x)
{
	long i;
	long str_beg, n_elem_str;
	double s;

	for(i=0; i<n; i++)
	{
		str_beg = ig[i];
		n_elem_str = ig[i+1] - str_beg;

		s = Scal(&ggl[str_beg], &x[i-n_elem_str], n_elem_str);

		if(fabs(di[i]) < this->zero)
		{
			printf("Error: Solve_L failed!!! Divizion by zero!!!\n");
			return 1;
		}

		x[i] = (pr[i] - s)/di[i];
	}
	return 0;
}
//-----------------------------------------------------------------
