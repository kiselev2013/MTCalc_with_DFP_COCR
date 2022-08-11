/**
 * GENERAL REMARKS
 *
 *  This code is freely available under the following conditions:
 *
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer.
 *
 *              MTCalc_with_DFP_COCR
 *  Header file for pardiso.cpp
 *
 *  Written by Prof. Marina G. Persova and Ph.D. Petr A. Domnikov 
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  mpersova@mail.ru
 *  Version 1.5 December 17, 2020
*/
#pragma once
#include "FormatConverter.h"

struct pardiso_solver
{
	FormatConverter f;

	// memmory for matrix
	MKL_INT64 *ia;
	MKL_INT64 *ja;
	double *a;

	// PARDISO parameters
	MKL_INT64 n;
	MKL_INT64 mtype;
	MKL_INT64 nrhs;
	void *pt[64];
	MKL_INT64 maxfct;
	MKL_INT64 mnum;
	MKL_INT64 msglvl;
	MKL_INT64 phase;
	MKL_INT64 *perm;
	MKL_INT64 iparm[64];
	MKL_INT64 info;

	pardiso_solver();
	~pardiso_solver();

	void factorize(int nb,int *ig,int *jg,double *ggl,double *ggu,double *di,int *idi,int *ijg,int nthreads);
	void solve_nrhs(int nrhs,double *pr,double *x);
	void stop_solver();
	void clear();
	void init();
	void init_rsf();
	void factorize_rsf(int nb,int *ig,int *jg,double *ggl,double *di,int nthreads);
};
