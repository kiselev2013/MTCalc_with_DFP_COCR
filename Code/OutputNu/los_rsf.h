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