#pragma once
#include "base_solver.h"

class RSF_solver : public Base_solver
{
/*
����� RSF_solver �������� ������� ������������ ���
�������� � ��������� � ����������� �������� �������
*/
public:
	RSF_solver();
	~RSF_solver();


	/////////////////////   �������� ������������ ��� ���   ///////////////////////////////

	int LLT(long *ig, long *jg, double *gg, double *di, double *sg, double *d, long n);

	int LLTd1(long *ig, long *jg, double *gg, double *di, double *sg, double *d, long n);

	int LU_sq(long *ig, long *jg, double *ggl, double *ggu, double *di,
		double *sl, double *su, double *d, long n);

	// 1 - � ������� L
	int LU(long *ig, long *jg, double *ggl, double *ggu, double *di,
		double *sl, double *su, double *d, long n);



	/////////////////////  ��������-��������� ���������  ////////////////////////////

	// ��������� ������������ ����������� ������� �� ������
	void mult_symmetr(long *ig, long *jg, double *gg, double *di, double *x, double *y, long n);

	// ��������� �������������� ����������� ������� �� ������
	void mult_mv(long *ig, long *jg, double *ggl, double *ggu, double *di, double *x, double *y, long n);
	void mult_mv_omp(long *ig, long *jg, double *ggl, double *ggu, double *di, double *x, double *y, long n, double *y_omp);

	// ��������� �� ����������������� ������� (�� ��������� �� 1)
	void mult_u_d(long *ig, long *jg, double *ggu, double *di, double *x, double *y, long n);

	// ��������� �� ���������������� ������� (�� ��������� �� 1)
	void mult_l_d(long *ig, long *jg, double *ggl, double *di, double *x, double *y, long n);


	/////////////////// ������� ���� � ������������ ��������� ///////////////////////

	// ������� ���� � ���������������� �������� (�� ��������� �� 1)
	void solve_l_d(long *ig, long *jg, double *ggl, double *di, double *f, double *x, long n);

	// ������� ���� � ����������������� �������� (�� ��������� �� 1)
	void solve_u_d(long *ig, long *jg, double *ggu, double *di, double *f, double *x, double *s, long n);


	//////////////////   ��������������� ��������   ///////////////////////////

	int Arnoldi(long n, long *ig, long *jg, double *di, double *ggl, double *ggu,
		double *v, double *h, double *w, long m, long *p);

	int Arnoldi_diag(long n, long *ig, long *jg, double *di, double *ggl, double *ggu,
		double *v, double *h, double *w, long m, long *p,
		double *d, double *help);

	int Arnoldi_lu_sq(long n, long *ig, long *jg, double *di, double *ggl, double *ggu,
		double *v, double *h, double *w, long m, long *p,
		double *d, double *help, double *sl, double *su);

	int Arnoldi_3diag(long n, long *ig, long *jg, double *di, double *ggl, double *ggu,
		double *v, double *h, double *w, long m, long *p,
		double *d, double *help, double *sl, double *su, long *ig_d, long *jg_d);





};