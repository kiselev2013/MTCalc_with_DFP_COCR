#include "StdAfx.h"
#include "rci.h"
extern ofstream logfile;
//------------------------------------------------------------------------
RCI::RCI(int n, int maxiter, double eps, double *x, double *pr, double **in, double **out)
{
	this->n = n;
	this->maxiter = maxiter;
	this->eps = eps;
	this->x = x;
	this->pr = pr;
	this->in = in;
	this->out = out;

	stage = 1;
	iter = 0;
	r_old = -1;
	eps_x0 = 1e-30;
	eps_zero = 1e-40;
	eps_gmres = 1e-20;
	residual = -1;
	residualRel = 1;
}
//------------------------------------------------------------------------
RCI::~RCI()
{
}
//------------------------------------------------------------------------
void RCI::PrintIterResidual()
{
	cout << iter << "    " << scientific << residualRel << endl << flush;
	logfile << iter << "    " << scientific << residualRel << endl << flush;
}
//------------------------------------------------------------------------
void RCI::DoStopTest(double *r, int *req)
{
	residual = bs.Norm_Euclid(r, n);
	residualRel = residual/r_old;

	if (residualRel < eps)
	{
		*req = REQ_OK;
		stage = 1;
	}
}
//------------------------------------------------------------------------
void RCI::DoStopTest(double r, int *req)
{
	residual = r;
	residualRel = residual/r_old;

	if (residualRel < eps)
	{
		*req = REQ_OK;
		stage = 1;
		PrintIterResidual();
	}
}
//------------------------------------------------------------------------
// проверяем, является ли уже начальное приближение решением
//------------------------------------------------------------------------
void RCI::DoX0Test(double *r, int *req)
{
	residual = bs.Norm_Euclid(r, n);

	if (residual < eps_x0)
	{
		logfile << "x0 is solution.\n";

		*req = REQ_X0_NULL;
		stage = 1;
	}
}	
//------------------------------------------------------------------------
