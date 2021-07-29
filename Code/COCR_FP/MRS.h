#pragma once
//------------------------------------------------------------------------
class MRS
{
private:
	int n;
	int nb;
	double *y;
	double *s;
	double *u;
	double w;
	double w_re;
	double w_im;
	double residual;
	bool isInit;
	double eps;

public:
	void Mrs(double *x, double *r);
	void MrsCmplx(double *x, double *r);

	double* GetResidual();
	void GetResidual(double *r);
	double GetResidualNorm();
	void GetSolution(double *a);
	double GetW();
	double GetWRe();
	double GetWIm();

	void CalcResidual();

	MRS(int n);
	~MRS();
};


