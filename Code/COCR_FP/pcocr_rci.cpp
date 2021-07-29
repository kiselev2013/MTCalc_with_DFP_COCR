#include "stdafx.h"
#include "pcocr_rci.h"
#include "ControlOMP.h"
extern ControlOMP omp;
extern ofstream logfile;
//------------------------------------------------------------------------
PCOCR_RCI::PCOCR_RCI(int n, int maxiter, double eps, double *x, double *pr, double **in, double **out) :
	RCI(n, maxiter, eps, x, pr, in, out)
{
	nb = n/2;

	r = new double[n];
	p = new double[n];
	s = new double[n];
	z = new double[n];
	a = new double[n];
	w = new double[n];
}
//------------------------------------------------------------------------
PCOCR_RCI::~PCOCR_RCI()
{
	if (r) {delete [] r; r=NULL;}
	if (p) {delete [] p; p=NULL;}
	if (s) {delete [] s; s=NULL;}
	if (z) {delete [] z; z=NULL;}
	if (a) {delete [] a; a=NULL;}
	if (w) {delete [] w; w=NULL;}
}
//------------------------------------------------------------------------
int PCOCR_RCI::Run()
{
	int i;

	switch(stage)
	{
	case 1:
		// вычисляем невязку системы (r)
		*in = x;
		*out = r;
		request = REQ_MULT_MV;
		stage = 2;
		break;
	case 2:
		#pragma omp parallel for private(i) num_threads(omp.GetNumberOfThreads())
		for(i=0; i<n; i++)
		{
			r[i] = pr[i] - r[i];
		}

		// норма невязки перед началом счёта
		// как будто считаем с нулевого начального приближения
		r_old = bs.Norm_Euclid(pr, n);
		//cout<<"NORM= "<<bs.Norm_Euclid(r, n)/bs.Norm_Euclid(pr, n)<<endl;

		// проверяем, является ли уже начальное приближение решением
		*in = r;
		request = REQ_X0_TEST;
		stage = 3;
		break;

	case 3:
		// невязка предобусловленной системы (s=M^{-1}*r)
		*in = r;
		*out = s;
		request = REQ_PRECOND;
		stage = 4;
		break;

	case 4:
		// p_{-1} = a_{-1} = 0
		#pragma omp parallel for private(i) num_threads(omp.GetNumberOfThreads())
		for (i=0; i<n; i++)
		{
			p[i] = a[i] = w[i] = 0;
		}

		// z=As
		*in = s;
		*out = z;
		request = REQ_MULT_MV;
		stage = 5;
		break;

	case 5:
MAIN_CYCLE:
		// главный цикл
		iter++;

		if(iter > maxiter)
		{
			request = REQ_MAXITER;
			stage = 1;
			break;
		}

		// p = s + betta*p
		bs.Cmplx_axpy(betta, s, p, p, nb);
		
		// a = z + betta*a 
		bs.Cmplx_axpy(betta, z, a, a, nb);

		// (z, s)
		zs = bs.ScalCmplx(z, s, nb);

		// w
		*in = a;
		*out = w;
		request = REQ_PRECOND;
		stage = 6;
		break;

	case 6:
		// (a, w)
		aw = bs.ScalCmplx(a, w, nb);
		if (fabs(real(aw)) < eps_zero  && fabs(imag(aw)) < eps_zero)
		{
			cout << "aw == 0\n";
			request = REQ_DIV_BY_ZERO;
			break;
		}

		// alpha
		alpha = zs/aw;

		// x = x + alpha*p
		bs.Cmplx_axpy(alpha, x, p, x, nb);

		// r = r - alpha*a
		bs.Cmplx_axpy(-alpha, r, a, r, nb);

		*in = r;
		request = REQ_STOP_TEST;
		stage = 7;
		break;

	case 7:
		// s = s - alpha*w
		bs.Cmplx_axpy(-alpha, s, w, s, nb);

		// z
		*in = s;
		*out = z;
		request = REQ_MULT_MV;
		stage = 8;
		break;

	case 8:
		//zs
		zs_1 = zs;

		zs = bs.ScalCmplx(z, s, nb);
		if (fabs(real(aw)) < eps_zero  && fabs(imag(aw)) < eps_zero)
		{
			cout << "zs == 0\n";
			request = REQ_DIV_BY_ZERO;
			break;
		}
		betta = zs/zs_1;	

		stage = 5;
		goto MAIN_CYCLE;
	}

	return request;
}
//------------------------------------------------------------------------
