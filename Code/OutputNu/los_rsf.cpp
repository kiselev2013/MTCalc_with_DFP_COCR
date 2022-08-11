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

#include "stdafx.h"
#include "los_rsf.h"
#include "ControlOMP.h"

extern ControlOMP omp;
extern ofstream logfile; 

//------------------------------------------------------------------------------------------
LOS_rsf_solver::LOS_rsf_solver()
{
	ShowProgress=1;
}
//------------------------------------------------------------------------------------------
LOS_rsf_solver::~LOS_rsf_solver()
{
}
//------------------------------------------------------------------------------------------
long LOS_rsf_solver::LOS_LU_sq(long n, long *ig, long *jg, double *di, double *ggl, double *ggu,
					   double *f, double *x, double eps, long maxiter,
					   double *d, double *sl, double *su,
					   double *r, double *z, double *p, double *qr, double *laqr, double *h)
{
	long i, iter;
	double alpha, beta;
	double r_old, r_new, r_nev;
	double p_scal;

	double *y_omp = new double[n*(omp.GetMaxThreads()-1)];

	logfile << "LOS_LU_sq...\n";

#if defined(__PROGRESS_INDICATOR_INCLUDED)
	CString sbuf;
	GUI::ProgressIndicator * pin=NULL;
	if(ShowProgress){
	pin=new GUI::ProgressIndicator("LOS_LU(sq)...", long(-log10(eps)*1000));
	}
#endif


	// initial approximation
	#pragma omp parallel shared(x) private(i) num_threads(omp.GetNumberOfThreads())
	{
		#pragma omp for
		for(i=0; i<n; i++)
		{
			x[i] = 0.0;
		}
	}

	// calculate the residual of the original system
	mult_mv_omp(ig, jg, ggl, ggu, di, x, h, n, y_omp);

	#pragma omp parallel shared(h, f) private(i) num_threads(omp.GetNumberOfThreads())
	{
		#pragma omp for
		for(i=0; i<n; i++)
		{
			h[i] = f[i] - h[i];
		}
	}

	// norm of the true residual  before the start of counting
	r_old = Norm_Euclid(h, n);

	// check if the initial approximation is already a solution
	if (r_old < 1e-30)
	{
		logfile << "x0 is solution.\n";
		iter = 0;
		r_nev = 0;
		goto label1;
	}

	// residual of the preconditioned system
	solve_l_d(ig, jg, sl, d, h, r, n);

	// z0
	solve_u_d(ig, jg, su, d, r, z, h, n);

	// p0
	mult_mv_omp(ig, jg, ggl, ggu, di, z, h, n, y_omp);
	solve_l_d(ig, jg, sl, d, h, p, n);

	/ iterative process
	for(iter=1; iter<=maxiter; iter++)
	{
		// alpha
		p_scal = Scal(p, p, n);
		alpha = Scal(p, r, n)/p_scal;

		// x, r
		#pragma omp parallel shared(x, r, z, p) private(i) num_threads(omp.GetNumberOfThreads())
		{
			#pragma omp for
			for (i=0; i<n; i++)
			{
				x[i] += alpha*z[i];
				r[i] -= alpha*p[i];
			}
		}

		// check for exit from the iterative process
		//solve_l_d(ig, jg, sl, d, r, h, n);
		mult_l_d(ig, jg, sl, d, r, h, n);

		r_new = Norm_Euclid(h, n);
		r_nev = r_new/r_old;

		logfile << iter << '\t' << scientific << r_nev << '\n';

#if defined(__PROGRESS_INDICATOR_INCLUDED)
		if(ShowProgress){
		sbuf.Format("LOS_LU(sq): eps = %e, eps current = %e, iter = %d", eps, r_nev, iter);
		pin->SetText(sbuf);	
		pin->SetPos(long(-log10(r_nev)*1000));
		}
#endif 

		if(r_nev < eps)
			goto label1;

		// beta
		solve_u_d(ig, jg, su, d, r, qr, h, n);
		mult_mv_omp(ig, jg, ggl, ggu, di, qr, h, n, y_omp);
		solve_l_d(ig, jg, sl, d, h, laqr, n);

		beta = -Scal(p, laqr, n)/p_scal;

		// z, p
		#pragma omp parallel shared(z, qr, p, laqr) private(i) num_threads(omp.GetNumberOfThreads())
		{
			#pragma omp for
			for (i=0; i<n; i++)
			{
				z[i] = qr[i] + beta*z[i];
				p[i] = laqr[i] + beta*p[i];
			}
		}
	}

label1:	cout << "LOS_LU_sq: iter="<<iter<<" residual="<<r_nev<<endl;
#if defined(__PROGRESS_INDICATOR_INCLUDED)
	if(ShowProgress){
	delete pin;
	}
#endif

	if (y_omp) {delete [] y_omp; y_omp=NULL;}

	return iter;
}
//------------------------------------------------------------------------------------------
//-- LOS with diagonal preconditioning and additional stopping criterion for a problem with a loop on vector elements 
//------------------------------------------------------------------------------------------
long LOS_rsf_solver::LOS_diag(long n, long *ig, long *jg, double *di, double *ggl, double *ggu,
					  double *f, double *x, double eps, long maxiter,
					  double *d, double *r, double *z, double *p,
					  double *qr, double *laqr, double *h, double *temp_x)
{
	long i, iter=0;
	double alpha, beta;
	double r_old, r_new, r_nev;
	double p_scal;
	double temp;

	// for early exit from the iterative process
	double r1, r2;
	long flag = 0;
	bool flag_x0=false;
	r1 = 1.0;
	double ssf; // right side vector norm

	logfile << "LOS_diag ...\n";

	// initial approximation
	for(i=0; i<n; i++)
		x[i] = 0.0;

	// diagonal preconditioner
	// here d[] are the roots of the modules of the diagonal elements to the minus one degree
	for(i=0; i<n; i++)
	{
		temp = fabs(di[i]);

		if(temp < 1e-30)
		{
			d[i] = 1.0; 
		}
		else
		{
			d[i] = 1.0/sqrt(temp);
		}
	}

	//  calculate the residual of the original system
	mult_mv(ig, jg, ggl, ggu, di, x, h, n);

	for(i=0; i<n; i++)
		h[i] = f[i] - h[i];

	// the norm of the true residual before the start of counting
	r_old = Norm_Euclid(h, n);

	// right-hand side norm
	ssf = this->Norm_Euclid(f, n);
	logfile << "|r|=" << scientific << r_old << "\t|f|=" << ssf << '\n';

	// check if the initial approximation is already a solution
	if (r_old < 1e-30)
	{
		logfile << "x0 is solution.\n";
		iter = 0;
		r_nev = r_old;
		goto label1;
	}

	if(r_old/ssf < eps)
	{
		logfile << "x0 is solution.\n";
		iter = 0;
		r_nev = r_old/ssf;
		goto label1;
	}

	if (r_old > ssf)
	{
		// this initial approximation is worse than zero
		// start from zero initial approximation
		logfile << "|r0|/|f|=" << scientific << r_old/ssf << " Starting with x0=0.\n";

		for (i=0; i<n; i++)
		{
			h[i] = f[i];
			x[i] = 0.0;
		}

		r_old = ssf;
	}

	// as if we are counting from zero initial approximatio
	r_old = ssf;

	// residual of the preconditioned system
	for(i=0; i<n; i++)
		r[i] = h[i]*d[i];

	// z0
	for(i=0; i<n; i++)
		z[i] = r[i]*d[i];

	// p0
	mult_mv(ig, jg, ggl, ggu, di, z, h, n);
	for(i=0; i<n; i++)
		p[i] = h[i]*d[i];

	//iterative process
	for(iter=1; iter<=maxiter; iter++)
	{
		// alpha
		p_scal = Scal(p, p, n);
		alpha = Scal(p, r, n)/p_scal;

		// x, r
		for (i=0; i<n; i++)
		{
			x[i] += alpha*z[i];
			r[i] -= alpha*p[i];
		}

		// checking for exit from the iterative process
		for(i=0; i<n; i++)
			h[i] = r[i]/d[i];

		r_new = Norm_Euclid(h, n);
		r_nev = r_new/r_old;

		logfile << iter << '\t' << scientific << r_nev << '\n';

		if(r_nev < eps)
			goto label1;

		// additional stop criterion
		if (r_nev < r1 && flag == 0)
		{
			for(i=0; i<n; i++)
				temp_x[i] = x[i];

			flag = 1;

			r2 = r_nev/sqrt(10.0);
		}
		else if (r_nev < r2 && flag == 1)
		{
			double norm_x, norm_xk_x;
			double razn;

			norm_x = Norm_Euclid(x, n);
			for(i=0; i<n; i++)
				temp_x[i] -= x[i];
			norm_xk_x = Norm_Euclid(temp_x, n);

			razn = norm_xk_x/norm_x;
			logfile << "razn=" << scientific << razn << '\n';

			for(i=0; i<n; i++)
				temp_x[i] = x[i];

			if (razn < 1e-6)
			{
				logfile << "Solution has not been changing -> exit.\n";
				break;
			}

			r2 /= sqrt(10.0);
		}

		// beta
		for(i=0; i<n; i++)
			qr[i] = r[i]*d[i];

		mult_mv(ig, jg, ggl, ggu, di, qr, h, n);
		for(i=0; i<n; i++)
			laqr[i] = h[i]*d[i];

		beta = -Scal(p, laqr, n)/p_scal;

		// z, p
		for (i=0; i<n; i++)
		{
			z[i] = qr[i] + beta*z[i];
			p[i] = laqr[i] + beta*p[i];
		}

	}

label1:	r[0] = r_nev;

	return iter;
}
//------------------------------------------------------------------------------------------