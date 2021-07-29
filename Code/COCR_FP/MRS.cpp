#include "StdAfx.h"
#include "MRS.h"
#include "base_solver.h"
#include "ControlOMP.h"
extern ControlOMP omp;

//------------------------------------------------------------------------
MRS::MRS(int n)
{
	this->n = n;
	nb = n/2;

	y = new double[n];
	s = new double[n];
	u = new double[n];

	isInit = false;

	eps = 1e-100;
}
//------------------------------------------------------------------------
MRS::~MRS(void)
{
	if(y) {delete [] y; y = NULL;}
	if(s) {delete [] s; s = NULL;}
	if(u) {delete [] u; u = NULL;}
}
//------------------------------------------------------------------------
double MRS::GetResidualNorm()
{
	return residual;
}
//------------------------------------------------------------------------
void MRS::GetResidual(double *r)
{
	memcpy_s(r, n*sizeof(double), s, n*sizeof(double));
}
//------------------------------------------------------------------------
double* MRS::GetResidual()
{
	return s;
}
//------------------------------------------------------------------------
double MRS::GetW()
{
	return w;
}
//------------------------------------------------------------------------
double MRS::GetWRe()
{
	return w_re;
}
//------------------------------------------------------------------------
double MRS::GetWIm()
{
	return w_im;
}
//------------------------------------------------------------------------
void MRS::GetSolution(double *a)
{
	memcpy_s(a, n*sizeof(double), y, n*sizeof(double));
}
//------------------------------------------------------------------------
void MRS::Mrs(double *x, double *r)
{
	int i;
	double su = 0.0;
	double uu = 0.0;

	if (!isInit)
	{
		memcpy_s(y, n*sizeof(double), x, n*sizeof(double));
		memcpy_s(s, n*sizeof(double), r, n*sizeof(double));

		CalcResidual();

		isInit = true;
	} 
	else
	{
		#pragma omp parallel shared(r) private(i) num_threads(omp.GetNumberOfThreads())
		{
			#pragma omp for reduction(+:su) reduction(+:uu)
			for (i=0; i<n; i++)
			{
				u[i] = r[i] - s[i];
				su += s[i]*u[i];
				uu += u[i]*u[i];
			}
		}

		if (fabs(uu) < eps)
		{
			w = 0;
		} 
		else
		{
			w = -su/uu;
		}

		if (w < 0)
		{
			w = 0;
		} 
		else if (w > 1) // проверить правда ли быстрее сходится
		{
			w = 1;

			memcpy_s(y, n*sizeof(double), x, n*sizeof(double));
			memcpy_s(s, n*sizeof(double), r, n*sizeof(double));

			CalcResidual();
		}
		else
		{
			double res = 0;

			#pragma omp parallel shared(x) private(i) num_threads(omp.GetNumberOfThreads())
			{
				#pragma omp for reduction(+:res)
				for (i=0; i<n; i++)
				{
					y[i] += w*(x[i] - y[i]);
					s[i] += w*u[i];
					res += s[i]*s[i];
				}
			}

			residual = sqrt(res);
		}
	}
}
//------------------------------------------------------------------------
void MRS::MrsCmplx(double *x, double *r)
{
	int i, i_re, i_im;
	double su_re = 0;
	double su_im = 0;
	double uu_re = 0;
	double uu_im = 0;
	std::complex<double> t;
	double t_re, t_im;

	if (!isInit)
	{
		memcpy_s(y, n*sizeof(double), x, n*sizeof(double));
		memcpy_s(s, n*sizeof(double), r, n*sizeof(double));

		CalcResidual();

		isInit = true;
	} 
	else
	{
		for (i=0; i<nb; i++)
		{
			i_re = i*2;
			i_im = i*2+1;

			u[i_re] = r[i_re] - s[i_re];
			u[i_im] = r[i_im] - s[i_im];

			su_re += s[i_re]*u[i_re] - s[i_im]*u[i_im];
			su_im += s[i_re]*u[i_im] + s[i_im]*u[i_re];

			uu_re += u[i_re]*u[i_re] - u[i_im]*u[i_im];
			uu_im += u[i_re]*u[i_im] + u[i_im]*u[i_re];
		}

		if (fabs(uu_re) < eps && fabs(uu_im) < eps)
		{
			w_re = w_im = 0;
		} 
		else
		{
			t = -std::complex<double>(su_re, su_im)/std::complex<double>(uu_re, uu_im);
			w_re = real(t); 
			w_im = imag(t); 
		}

		if (fabs(w_re) < eps && fabs(w_im) < eps)
		{
			w_re = w_im = 0;
		} 
		else if (abs(std::complex<double>(w_re, w_im)) > 1) // проверить правда ли быстрее сходится
		{
			w_re = 1;
			w_im = 0;

			memcpy_s(y, n*sizeof(double), x, n*sizeof(double));
			memcpy_s(s, n*sizeof(double), r, n*sizeof(double));

			CalcResidual();
		}
		else
		{
			residual = 0;

			for (i=0; i<nb; i++)
			{
				i_re = i*2;
				i_im = i*2+1;

				t_re = x[i_re] - y[i_re];
				t_im = x[i_im] - y[i_im];

				y[i_re] += w_re*t_re - w_im*t_im;
				y[i_im] += w_re*t_im + w_im*t_re;
				
				s[i_re] += w_re*u[i_re] - w_im*u[i_im];
				s[i_im] += w_re*u[i_im] + w_im*u[i_re];

				residual += s[i_re]*s[i_re] + s[i_im]*s[i_im];
			}

			residual = sqrt(residual);
		}
	}
}
//------------------------------------------------------------------------
void MRS::CalcResidual()
{
	Base_solver a;
	residual = a.Norm_Euclid(s, n);
}
//------------------------------------------------------------------------
