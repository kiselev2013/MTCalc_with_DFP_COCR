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
                                                                                                                                                                                                                                                                                            
#include "stdafx.h"
#include "rsf_solver.h"

#include "ControlOMP.h"
extern ControlOMP omp;

//------------------------------------------------------------------------------------------
RSF_solver::RSF_solver()
{
}
//------------------------------------------------------------------------------------------
RSF_solver::~RSF_solver()
{
}
//------------------------------------------------------------------------------------------
void RSF_solver::mult_symmetr(long *ig, long *jg, double *gg, double *di, double* x, double *y, long n)
{
	mult_mv(ig, jg, gg, gg, di, x, y, n);
}
//-----------------------------------------------------------------------
void RSF_solver::mult_mv(long *ig, long *jg, double *ggl, double *ggu,
						 double *di, double *x, double *y, long n)
{
	long i, j, k;

	for (i=0; i<n; i++)
	{
		y[i] = 0.0;
	}

	for(i=0; i<n; i++)
	{
		y[i] = di[i]*x[i];
		for(j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			y[i] += ggl[j]*x[k];
			y[k] += ggu[j]*x[i];
		}
	}
}
//------------------------------------------------------------------------
void RSF_solver::mult_mv_omp(long *ig, long *jg, double *ggl, double *ggu, double *di, double *x, double *y, long n, double *y_omp)
{
	long i, j, k;
	int rank, adr;
	int nThreads = omp.GetNumberOfThreads();

	#pragma omp parallel shared(ig, jg, ggl, ggu, di, x, y, y_omp) private(i, j, k, rank, adr) num_threads(nThreads)
	{
		#pragma omp for nowait
		for (i=0; i<n; i++)
		{
			y[i] = 0.0;
		}

		#pragma omp for
		for (i=0; i<n*(nThreads-1); i++)
		{
			y_omp[i] = 0.0;
		}

		#pragma omp for
		for(i=0; i<n; i++)
		{
			rank = 0;

			if (rank == 0)
			{
				y[i] = di[i]*x[i];
				for(j=ig[i]; j<=ig[i+1]-1; j++)
				{
					k = jg[j];
					y[i] += ggl[j]*x[k];
					y[k] += ggu[j]*x[i];
				}
			}
			else
			{
				adr = (rank - 1)*n;
				y_omp[adr + i] = di[i]*x[i];

				for(j=ig[i]; j<=ig[i+1]-1; j++)
				{
					k = jg[j];
					y_omp[adr + i] += ggl[j]*x[k];
					y_omp[adr + k] += ggu[j]*x[i];
				}
			}
		}//i
		//reduction y_omp

		#pragma omp for
		for (i=0; i<n; i++)
		{
			for (j=0; j<nThreads-1; j++)
				y[i] += y_omp[j*n + i];
		}	
	}//parallel
}
//------------------------------------------------------------------------------------------
int RSF_solver::LLT(long *ig, long *jg, double *gg, double *di, double *sg, double *d, long n)
{
	long i, l, k, k1;
	double s, temp;

	printf("LLT-decomposition... ");
	for(i=0; i<n; i++)
		d[i] = 0.0;
	for(i=0; i<ig[n]; i++)
		sg[i] = 0.0;

	for(i=0; i<n; i++)
	{
		for(l=ig[i]; l<=ig[i+1]-1; l++)
		{
			s = 0.0;
			for(k=ig[i]; k<=l-1; k++)
				for(k1=ig[jg[l]]; k1<=ig[jg[l]+1]-1; k1++)
					if(jg[k1] == jg[k])
					{
						s += sg[k1]*sg[k];
						break;
					}
					sg[l] = (gg[l] - s)/d[jg[l]];
		}
		s = 0.0;
		for(k=ig[i]; k<=ig[i+1]-1; k++)
			s += sg[k]*sg[k];
		temp = di[i] - s;
		if(temp >= 0.0)
		{
			d[i] = sqrt(temp);
		}
		else
		{
			printf("LLT failed i=%ld n=%ld.\n", i, n);
			return 1;
		}
	}
	printf("done.\n");
	return 0;
}
//------------------------------------------------------------------------------------------
int RSF_solver::LLTd1(long *ig, long *jg, double *gg, double *di, double *sg, double *d, long n)
{
	long i, l, k, k1;
	double s, temp;
	std::vector<double> d1;

	d1.resize(n);

	printf("LLT-decomposition... ");
	for(i=0; i<n; i++)
		d[i] = 0.0;
	for(i=0; i<ig[n]; i++)
		sg[i] = 0.0;

	for(i=0; i<n; i++)
	{
		for(l=ig[i]; l<=ig[i+1]-1; l++)
		{
			s = 0.0;
			for(k=ig[i]; k<=l-1; k++)
			{
				for(k1=ig[jg[l]]; k1<=ig[jg[l]+1]-1; k1++)
					if(jg[k1] == jg[k])
					{
						s += sg[k1]*d1[jg[k]]*sg[k];
						break;
					}
			}

			if(fabs(d[jg[l]])<1e-30)
			{
				printf("LLT failed i=%ld n=%ld.\n", i, n);
				return 1;
			}
			sg[l] = (gg[l] - s)/(d[jg[l]]*d1[jg[l]]);
		}

		s = 0.0;
		for(k=ig[i]; k<=ig[i+1]-1; k++)
			s += sg[k]*d1[jg[k]]*sg[k];
		temp = di[i] - s;


		if(temp < 0.0)
		{
			d1[i] = -1;
		}
		else
		{
			d1[i] = 1;
		}

		d[i] = sqrt(fabs(temp));
	}
	printf("done.\n");
	return 0;
}
//----------------------------------------------------------
int RSF_solver:: LU_sq(long *ig, long *jg, double *ggl, double *ggu, double *di,
						double *sl, double *su, double *d, long n)
{
	printf("LU(sq)-decomposition... ");
	long i, l, k, k1;
	double s, temp;

	for(i=0; i<n; i++)
		d[i] = 0.0;

	for(i=0; i<ig[n]; i++)
		sl[i] = su[i] = 0.0;

	for(i=0; i<n; i++)
	{
		for(l=ig[i]; l<=ig[i+1]-1; l++)
		{
			s = 0.0;
			for(k=ig[i]; k<=l-1; k++)
				for(k1=ig[jg[l]]; k1<=ig[jg[l]+1]-1; k1++)
					if(jg[k1] == jg[k])
					{
						s += sl[k1]*su[k];
						break;
					}
					if(fabs(d[jg[l]])>1e-30)
					{
						su[l] = (ggu[l] - s)/d[jg[l]];
					}
					else
					{
						printf("LU(sq) failed!!! Division by zero!!!\n");
						return 1;
					}
		}

		for(l=ig[i]; l<=ig[i+1]-1; l++)
		{
			s = 0.0;
			for(k=ig[i]; k<=l-1; k++)
				for(k1=ig[jg[l]]; k1<=ig[jg[l]+1]-1; k1++)
					if(jg[k1] == jg[k])
					{
						s += su[k1]*sl[k];
						break;
					}
					if(fabs(d[jg[l]])>1e-30)
					{
						sl[l] = (ggl[l] - s)/d[jg[l]];
					}
					else
					{
						printf("LU(sq) failed!!! Division by zero!!!\n");
						return 1;
					}
		}

		s = 0.0;
		for(k=ig[i]; k<=ig[i+1]-1; k++)
			s += sl[k]*su[k];

		temp = di[i] - s;
		if(temp >= 0.0)
		{
			d[i] = sqrt(temp);
		}
		else
		{
			printf("LU_sq failed!!! sqrt!!!\n");
			return 1;
		}
	}

	printf("done.\n");
	return 0;
}
//---------------------------------------------------------------------
int RSF_solver:: LU(long *ig, long *jg, double *ggl, double *ggu, double *di,
				   double *sl, double *su, double *d, long n)
{
	printf("LU-decomposition... ");
	int i, l, k, k1;
	double s, temp;

	for(i=0; i<n; i++)
		d[i] = 0.0;

	for(i=0; i<ig[n]; i++)
		sl[i] = su[i] = 0.0;

	for(i=0; i<n; i++)
	{
		for(l=ig[i]; l<=ig[i+1]-1; l++)
		{
			s = 0.0;
			for(k=ig[i]; k<=l-1; k++)
				for(k1=ig[jg[l]]; k1<=ig[jg[l]+1]-1; k1++)
					if(jg[k1] == jg[k])
					{
						s += sl[k1]*su[k];
						break;
					}
					su[l] = ggu[l] - s;
		}

		for(l=ig[i]; l<=ig[i+1]-1; l++)
		{
			s = 0.0;
			for(k=ig[i]; k<=l-1; k++)
				for(k1=ig[jg[l]]; k1<=ig[jg[l]+1]-1; k1++)
					if(jg[k1] == jg[k])
					{
						s += su[k1]*sl[k];
						break;
					}
					if(fabs(d[jg[l]])>1e-30)
					{
						sl[l] = (ggl[l] - s)/d[jg[l]];
					}
					else
					{
						printf("LU failed!!! Division by zero!!!\n");
						return 1;
					}
		}

		s = 0.0;
		for(k=ig[i]; k<=ig[i+1]-1; k++)
			s += sl[k]*su[k];

		d[i] = di[i] - s;
	}

	printf("done.\n");
	return 0;
}
//---------------------------------------------------------------------
void RSF_solver::solve_l_d(long *ig, long *jg, double *ggl, double *di, double *f, double *x, long n)
{
	long i, j, k;
	double s;

	for(i=0; i<n; i++)
	{
		s = 0.0;
		for(j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			s += x[k]*ggl[j];
		}
		x[i] = (f[i] - s)/di[i];
	}
}
//---------------------------------------------------------- 
void RSF_solver::solve_u_d(long *ig, long *jg, double *ggu, double *di, double *f, double *x, double *s, long n)
{
	long i, j, k;

	for(i=0; i<n; i++)
		s[i] = 0.0;

	for(i=n-1; i>=0; i--)
	{
		x[i] = (f[i] - s[i])/di[i];
		for(j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			s[k] += x[i]*ggu[j];
		}
	}
}
//----------------------------------------------------------
//------------------------------------------------------------
//------------ Arnoldi orthogonalization ----------------------
//------------------------------------------------------------

//------ no precondition for Compressed Row Storage format ---------------
int RSF_solver::Arnoldi(long n, long *ig, long *jg, double *di, double *ggl, double *ggu,
						 double *v, double *h, double *w, long m, long *p)
{
	// p - the number of vectors that could be constructed (p<=m)
	// h(m+1,m); v(m+1,n)

	long i, j, k;
	long ind;
	double h_ij, h_j1_j, t;

	*p = m;

	for(j=0; j<m; j++)
	{
		mult_mv(ig, jg, ggl, ggu, di, &v[j*n], w, n);

		for(i=0; i<=j; i++)
		{
			ind = i*n;
			h_ij = h[i*m+j] = Scal(w, &v[ind], n);
			for(k=0; k<n; k++)
				w[k] -= h_ij*v[ind+k];
		}

		t = Scal(w,w,n);
		if(fabs(t) < 1e-10)
		{
			*p = j+1;
			return 1;
		}

		h_j1_j = h[(j+1)*m + j] = sqrt(t);	

		ind = (j+1)*n;
		for(k=0; k<n; k++)
			v[ind+k] = w[k]/h_j1_j;
	}

	return 0;
}
//---------------------------------------------------------------------------------------
int RSF_solver::Arnoldi_diag(long n, long *ig, long *jg, double *di, double *ggl, double *ggu,
							  double *v, double *h, double *w, long m, long *p,
							  double *d, double *help)
{
	// p - число векторов, к-рые удалось построить (p<=m)
	// h(m+1,m); v(m+1,n)

	long i, j, k;
	long ind;
	double h_ij, h_j1_j, t;

	*p = m;

	for(j=0; j<m; j++)
	{
		ind = j*n;
		for(i=0; i<n; i++)
			w[i] = v[ind+i]/d[i];

		mult_mv(ig, jg, ggl, ggu, di, w, help, n);

		for(i=0; i<n; i++)
			w[i] = help[i]/d[i];

		for(i=0; i<=j; i++)
		{
			ind = i*n;
			h_ij = h[i*m+j] = Scal(w, &v[ind], n);
			for(k=0; k<n; k++)
				w[k] -= h_ij*v[ind+k];
		}

		t = Scal(w,w,n);
		if(fabs(t) < 1e-10)
		{
			*p = j+1;
			return 1;
		}

		h_j1_j = h[(j+1)*m + j] = sqrt(t);	

		ind = (j+1)*n;
		for(k=0; k<n; k++)
			v[ind+k] = w[k]/h_j1_j;
	}

	return 0;
}
//---------------------------------------------------------------------------------------
int RSF_solver::Arnoldi_lu_sq(long n, long *ig, long *jg, double *di, double *ggl, double *ggu,
							   double *v, double *h, double *w, long m, long *p,
							   double *d, double *help, double *sl, double *su)
{
	// p - the number of vectors that could be constructed (p<=m)
	// h(m+1,m); v(m+1,n)

	long i, j, k;
	long ind;
	double h_ij, h_j1_j, t;

	*p = m;

	for(j=0; j<m; j++)
	{
		solve_u_d(ig, jg, su, d, &v[j*n], w, help, n);
		mult_mv(ig, jg, ggl, ggu, di, w, help, n);
		solve_l_d(ig, jg, sl, d, help, w, n);

		for(i=0; i<=j; i++)
		{
			ind = i*n;
			h_ij = h[i*m+j] = Scal(w, &v[ind], n);
			for(k=0; k<n; k++)
				w[k] -= h_ij*v[ind+k];
		}

		t = Scal(w,w,n);
		if(fabs(t) < 1e-20)
		{
			*p = j+1;
			return 1;
		}

		h_j1_j = h[(j+1)*m + j] = sqrt(t);	

		ind = (j+1)*n;
		for(k=0; k<n; k++)
			v[ind+k] = w[k]/h_j1_j;
	}

	return 0;
}
//---------------------------------------------------------------------------------------
int RSF_solver::Arnoldi_3diag(long n, long *ig, long *jg, double *di, double *ggl, double *ggu,
							   double *v, double *h, double *w, long m, long *p,
							   double *d, double *help, double *sl, double *su,
							   long *ig_d, long *jg_d)
{
	// p - the number of vectors that could be constructed (p<=m)
	// h(m+1,m); v(m+1,n)

	long i, j, k;
	long ind;
	double h_ij, h_j1_j, t;

	*p = m;

	for(j=0; j<m; j++)
	{
		solve_u_d(ig_d, jg_d, su, d, &v[j*n], w, help, n);
		mult_mv(ig, jg, ggl, ggu, di, w, help, n);
		solve_l_d(ig_d, jg_d, sl, d, help, w, n);

		for(i=0; i<=j; i++)
		{
			ind = i*n;
			h_ij = h[i*m+j] = Scal(w, &v[ind], n);
			for(k=0; k<n; k++)
				w[k] -= h_ij*v[ind+k];
		}

		t = Scal(w,w,n);
		if(fabs(t) < 1e-20)
		{
			*p = j+1;
			return 1;
		}

		h_j1_j = h[(j+1)*m + j] = sqrt(t);	

		ind = (j+1)*n;
		for(k=0; k<n; k++)
			v[ind+k] = w[k]/h_j1_j;
	}

	return 0;
}
//------------------------------------------------------------
void RSF_solver::mult_u_d(long *ig, long *jg, double *ggu, double *di, double *x, double *y, long n)
{
	long i, j, k;

	for(i=0; i<n; i++)
	{
		y[i] = di[i]*x[i];
		for(j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			y[k] += ggu[j]*x[i];
		}
	}	
}
//------------------------------------------------------------
void RSF_solver::mult_l_d(long *ig, long *jg, double *ggl, double *di, double *x, double *y, long n)
{
	long i, j;

	#pragma omp parallel shared(ig, jg, ggl, di, x, y) private(i, j) num_threads(omp.GetNumberOfThreads())
	{
		for(i=0; i<n; i++)
		{
			y[i] = di[i]*x[i];
			for(j=ig[i]; j<=ig[i+1]-1; j++)
				y[i] += ggl[j]*x[jg[j]];
		}
	}
}
//------------------------------------------------------------
