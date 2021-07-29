#include "stdafx.h"
#include "pardiso.h"

extern ofstream logfile;

pardiso_solver::pardiso_solver()
{
	ia=NULL;
	ja=NULL;
	a=NULL;
	perm=NULL;

	n = 0;

	init();
}

pardiso_solver::~pardiso_solver()
{
	clear();
}

void pardiso_solver::init()
{
	int i;

	mtype = 6; //2;
	nrhs = 1;
	maxfct = 1;
	mnum = 1;
	msglvl = 1;
	phase = 13;

	for(i=0;i<64;i++){pt[i]=0;}
	for(i=0;i<64;i++){iparm[i]=0;}

	msglvl=0;
}

void pardiso_solver::init_rsf()
{
	int i;

	mtype = 2;
	nrhs = 1;
	maxfct = 1;
	mnum = 1;
	msglvl = 1;
	phase = 13;

	for(i=0;i<64;i++){pt[i]=0;}
	for(i=0;i<64;i++){iparm[i]=0;}

	msglvl=0;
}

void pardiso_solver::clear()
{
	if(a){delete [] a;  a=NULL;}
	if(ia){delete [] ia; ia=NULL;}
	if(ja){delete [] ja; ja=NULL;}
	if(perm){delete [] perm; perm=NULL;}
}

void pardiso_solver::factorize(int nb,int *ig,int *jg,double *ggl,double *ggu,double *di,int *idi,int *ijg,int nthreads)
{
	int i;
	int ig_n_1=ig[nb];
	int sz_iptr=0;
	int sz_jptr=0;

	init();

	if (n == 0)
	{
		clear();
	}

	mkl_set_num_threads(nthreads);

	// ��������������� � CSR
	f.From2x2ToCSR_Complex_1_Sym(nb, ig, idi, ijg, &sz_iptr, &sz_jptr);

	if (ia==NULL) ia = new MKL_INT64[sz_iptr];
	if (ja==NULL) ja = new MKL_INT64[sz_jptr];
	if (a==NULL)   a = new double[sz_jptr*2];

	f.From2x2ToCSRComplex_2_Sym	(nb, ig, jg, idi, ijg, di, ggl, ia, ja, a);

	for(i=0;i<sz_iptr;i++){ia[i]++;}
	for(i=0;i<sz_jptr;i++){ja[i]++;}

	if (perm==NULL)	perm = new MKL_INT64[nb];


	phase = 12;
	n = nb;

	pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, NULL, NULL, &info);

	cout<<"info= "<<info<<'\n';
}

void pardiso_solver::solve_nrhs(int _nrhs,double *pr,double *x)
{
	int i;
	phase = 33;
	nrhs = _nrhs;
	pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, pr, x, &info);
}

void pardiso_solver::stop_solver()
{
	phase = -1;
	pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, NULL, NULL, &info);
}

void pardiso_solver::factorize_rsf(int nb,int *ig,int *jg,double *ggl,double *di,int nthreads)
{
	int i;
	int ig_n_1=0;
	int sz_iptr=0;
	int sz_jptr=0;

	init_rsf();

	if(n!=nb)
	{
		clear();
	}

	mkl_set_num_threads(nthreads);

	// ��������������� � CSR
	f.FromRSFToCSR_Real_1_Sym(nb, ig, &sz_iptr, &sz_jptr);

	if(n!=nb)
	{
		ia = new MKL_INT64[sz_iptr];
		ja = new MKL_INT64[sz_jptr];
		a = new double[sz_jptr];
	}

	f.FromRSFToCSR_Real_2_Sym(nb, ig, jg, di, ggl, ia, ja, a);

	for(i=0;i<sz_iptr;i++){ia[i]++;}
	for(i=0;i<sz_jptr;i++){ja[i]++;}

	if(n!=nb)
	{
		perm = new MKL_INT64[nb];
	}

	phase = 12;
	n = nb;

	pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, NULL, NULL, &info);
}
