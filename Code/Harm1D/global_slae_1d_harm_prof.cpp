#include "stdafx.h"
#include "global_slae_1d_harm_prof.h"
#include "local_matrix_1d.h"
//-----------------------------------------------------------
Global_slae_1d_harm_prof::Global_slae_1d_harm_prof(long n, long n_elem,
												   long *ig, long *nvkat,
												   double *mu, double *sigma,
												   double omega, double *xyz,
												   long ay_hy)
{
	this->n = n;
	this->n_elem = n_elem;
	this->ig = ig;
	this->nvkat = nvkat;
	this->mu = mu;
	this->sigma = sigma;
	this->omega = omega;
	this->xyz = xyz;

    this->ig_n_1 = ig[n];

	ggl = NULL;
	ggu = NULL;
	di = NULL;
	pr = NULL;

	ggl = new double[ig_n_1];
	if(ggl == 0)
		Memory_allocation_error("ggl", "Global_slae_1d_harm_prof::Global_slae_1d_harm_prof");

	ggu = new double[ig_n_1];
	if(ggu == 0)
		Memory_allocation_error("ggu", "Global_slae_1d_harm_prof::Global_slae_1d_harm_prof");

	di = new double[n];
	if(di == 0)
		Memory_allocation_error("di", "Global_slae_1d_harm_prof::Global_slae_1d_harm_prof");

	pr = new double[n];
	if(pr == 0)
		Memory_allocation_error("pr", "Global_slae_1d_harm_prof::Global_slae_1d_harm_prof");

	this->ay_hy = ay_hy;
}
//-----------------------------------------------------------
Global_slae_1d_harm_prof::~Global_slae_1d_harm_prof()
{
	if(di) {delete [] di; di=NULL;}
	if(pr) {delete [] pr; pr=NULL;}
	if(ggl) {delete [] ggl; ggl=NULL;}
	if(ggu) {delete [] ggu; ggu=NULL;}
}
//-----------------------------------------------------------
void Global_slae_1d_harm_prof::Assembling_for_1d_harm_problem()
{
	long i, j, k;
	double f_re[2] = {0.0, 0.0};
	double f_im[2] = {0.0, 0.0};
	double h;
	long str, col;
	long el_str; // количество ненулевых эл-тов в текущей строке
	long el_0;   // количество нулевых эл-тов в текущей строке
	long adr;    // адрес элемента (str,col) в массиве ggl
	long displacement;

	// заполнение нулями
	for(i=0; i<n; i++)
		di[i] = pr[i] = 0.0;

	for(i=0; i<ig_n_1; i++)
		ggl[i] = ggu[i] = 0.0;

	// цикл по элементам
	for(i=0; i<n_elem; i++)
	{
		// вычисляем элементы локальной матрицы
		h = xyz[i+1] - xyz[i];
		Local_matrix_1d L(h, mu[nvkat[i]], sigma[nvkat[i]], omega, f_re, f_im, ay_hy);
		L.Compute_local_matrix_harm();

		// заносим элементы локальной матрицы в глобальную
		for(j=0; j<4; j++)
		{
			str = i*2 + j;
			di[str] += L.a[j][j];
			el_str = ig[str+1] - ig[str];
			el_0 = str - el_str;
			displacement = ig[str] - el_0;
			for(k=0; k<4; k++)
			{
				col = i*2 + k;
				if(col < str)
				{
					adr = displacement + col;
					ggl[adr] += L.a[j][k];
					ggu[adr] += L.a[k][j];
				}
			}
		}
	}
	
	if(ay_hy == 0)
	{
		Set_boundary_conditions_for_Ay();
	}
	else
	{
		Set_boundary_conditions_for_Hy();
	}
}
//-----------------------------------------------------------
void Global_slae_1d_harm_prof::Set_boundary_conditions_for_Ay()
{
	// бак
	di[0] = di[1] = 1.0;
	pr[0] = pr[1] = 0.0;
	ggl[0] = 0.0;	
	for(long i=0; i<5; i++)
		ggu[i] = 0.0;

	// второе краевое на поверхности
	pr[n-2] += 1.0;
}
//-----------------------------------------------------------
void Global_slae_1d_harm_prof::Set_boundary_conditions_for_Hy()
{
	long i, k;
	// бак
	di[0] = di[1] = 1.0;
	pr[0] = pr[1] = 0.0;
	ggl[0] = 0.0;	
	for(i=0; i<5; i++)
		ggu[i] = 0.0;

	// первое краевое на поверхности
	pr[n-2] = 1.0;
	pr[n-1] = 0.0;
	di[n-1] = di[n-2] = 1.0;

	for(i=0; i<5; i++)
	{
		k = ig_n_1-1-i;
		ggl[k] = 0.0;
	}

	ggu[ig_n_1-1] = 0.0;
}
//-----------------------------------------------------------