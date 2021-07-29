#pragma once
#include "rci.h"
//------------------------------------------------------------------------
class PCOCR_RCI : public RCI
{
private:
	int nb;
	std::complex<double> alpha, betta;
	std::complex<double> zs, aw, zs_1;
	double *r; // невязка исходной системы
	double *p; // направление поиска
	double *s; // невязка предобусловленной системы
	double *z;  // z=A*s
	double *a;  // a=A*p
	double *w;  // w=M^{-1}*a

public:
	PCOCR_RCI(int n, int maxiter, double eps, double *x, double *pr, double **in, double **out);
	~PCOCR_RCI();

	int Run();
};
