#include "stdafx.h"
#include "local_matrix_1d.h"
//-----------------------------------------------------------
Local_matrix_1d::Local_matrix_1d(double h, double mu, double sigma, double omega,
								 double *f_re, double *f_im, long alpha)
{
	this->h = h;
	this->mu = mu;
	this->sigma = sigma;
	this->omega = omega;

	for(long i=0; i<2; i++)
	{
		this->f_re[i] = f_re[i];
		this->f_im[i] = f_im[i];
	}

	this->alpha = alpha;
}
//-----------------------------------------------------------
Local_matrix_1d::~Local_matrix_1d()
{
}
//-----------------------------------------------------------
void Local_matrix_1d::Compute_local_matrix_b()
{
	double t;

	if(alpha==0)
	{
		t = 1.0/(mu*h);
	}
	else
	{
		t = 1.0/(sigma*h);
	}

	b[0][0] = b[1][1] = t;
    b[0][1] = b[1][0] = -t;
}
//-----------------------------------------------------------
void Local_matrix_1d::Compute_local_matrix_c()
{
	double t;

	if(alpha==0)
	{
		t = sigma*omega*h/6.0;
	}
	else
	{
		t = mu*omega*h/6.0;
	}

	c[0][0] = c[1][1] = 2.0*t;
    c[0][1] = c[1][0] = t;
}
//-----------------------------------------------------------
void Local_matrix_1d::Compute_local_matrix_harm()
{
	this->Compute_local_matrix_b();
	this->Compute_local_matrix_c();

	a[0][0] = b[0][0];
	a[0][1] = -c[0][0];
	a[0][2] = b[0][1];
	a[0][3] = -c[0][1];

	a[1][0] = c[0][0];
	a[1][1] = b[0][0];
	a[1][2] = c[0][1];
	a[1][3] = b[0][1];

	a[2][0] = b[1][0];
	a[2][1] = -c[1][0];
	a[2][2] = b[1][1];
	a[2][3] = -c[1][1];

	a[3][0] = c[1][0];
	a[3][1] = b[1][0];
	a[3][2] = c[1][1];
	a[3][3] = b[1][1];
}
//-----------------------------------------------------------
void Local_matrix_1d::Compute_local_vector(double *f, double *g)
{
	double t = h/6.0;

	g[0] = (2.0*f[0] + f[1])*t;
	g[1] = (f[0] + 2.0*f[1])*t;
}
//-----------------------------------------------------------
void Local_matrix_1d::Compute_local_vector_harm()
{
	Compute_local_vector(f_re, g_re);
	Compute_local_vector(f_im, g_im);

	g[0] = g_re[0];
	g[1] = g_im[0];
	g[2] = g_re[1];
	g[3] = g_im[1];
}
//-----------------------------------------------------------
void Local_matrix_1d::Compute_local_matrix_and_vector_harm()
{
	this->Compute_local_matrix_harm();
	this->Compute_local_vector_harm();
}
//----------------------------------------------------------------
