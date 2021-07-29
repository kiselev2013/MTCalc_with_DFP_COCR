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


	// начальное приближение
	#pragma omp parallel shared(x) private(i) num_threads(omp.GetNumberOfThreads())
	{
		#pragma omp for
		for(i=0; i<n; i++)
		{
			x[i] = 0.0;
		}
	}

	// вычисляем невязку исходной системы
	mult_mv_omp(ig, jg, ggl, ggu, di, x, h, n, y_omp);

	#pragma omp parallel shared(h, f) private(i) num_threads(omp.GetNumberOfThreads())
	{
		#pragma omp for
		for(i=0; i<n; i++)
		{
			h[i] = f[i] - h[i];
		}
	}

	// норма истинной невязки перед началом счёта
	r_old = Norm_Euclid(h, n);

	// проверяем, является ли уже начальное приближение решением
	if (r_old < 1e-30)
	{
		logfile << "x0 is solution.\n";
		iter = 0;
		r_nev = 0;
		goto label1;
	}

	// невязка предобусловленной системы
	solve_l_d(ig, jg, sl, d, h, r, n);

	// z0
	solve_u_d(ig, jg, su, d, r, z, h, n);

	// p0
	mult_mv_omp(ig, jg, ggl, ggu, di, z, h, n, y_omp);
	solve_l_d(ig, jg, sl, d, h, p, n);

	//итерационный процесс
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

		// проверка на выход из итерационного процесса
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
//-- ЛОС с диагональным предобусловливанием и дополнительным критерием останова
//-- для задачи с петлёй на векторных элементах
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

	// для досрочного выхода из итерационного процесса
	double r1, r2;
	long flag = 0;
	bool flag_x0=false;
	r1 = 1.0;
	double ssf; // норма вектора правой части

	logfile << "LOS_diag ...\n";

	// начальное приближение
	for(i=0; i<n; i++)
		x[i] = 0.0;

	// диагональный предобусловливатель
	// здесь d[] - корни из модулей диагональных эл-тов в минус первой
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

	// вычисляем невязку исходной системы
	mult_mv(ig, jg, ggl, ggu, di, x, h, n);

	for(i=0; i<n; i++)
		h[i] = f[i] - h[i];

	// норма истинной невязки перед началом счёта
	r_old = Norm_Euclid(h, n);

	// норма правой части
	ssf = this->Norm_Euclid(f, n);
	logfile << "|r|=" << scientific << r_old << "\t|f|=" << ssf << '\n';

	// проверяем, является ли уже начальное приближение решением
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
		// данное начальное приближение хуже нулевого
		// начинаем с нулевого начального приближения
		logfile << "|r0|/|f|=" << scientific << r_old/ssf << " Starting with x0=0.\n";

		for (i=0; i<n; i++)
		{
			h[i] = f[i];
			x[i] = 0.0;
		}

		r_old = ssf;
	}

	// как будто мы считаем с нулевого начального приближения
	r_old = ssf;

	// невязка предобусловленной системы
	for(i=0; i<n; i++)
		r[i] = h[i]*d[i];

	// z0
	for(i=0; i<n; i++)
		z[i] = r[i]*d[i];

	// p0
	mult_mv(ig, jg, ggl, ggu, di, z, h, n);
	for(i=0; i<n; i++)
		p[i] = h[i]*d[i];

	//итерационный процесс
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

		// проверка на выход из итерационного процесса
		for(i=0; i<n; i++)
			h[i] = r[i]/d[i];

		r_new = Norm_Euclid(h, n);
		r_nev = r_new/r_old;

		logfile << iter << '\t' << scientific << r_nev << '\n';

		if(r_nev < eps)
			goto label1;

		// дополнительный критерий останова
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