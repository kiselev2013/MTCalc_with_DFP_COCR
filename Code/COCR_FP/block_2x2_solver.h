#pragma once
#include "base_solver.h"
//------------------------------------------------------------------------

class Block_2x2_solver 
{
public:
	Block_2x2_solver();
	~Block_2x2_solver();

	// ��������� ������ �����
	inline void Mult_block_2x2(double *a, int size, double *x, double *y);
	inline void Mult_MV_block_2x2_transp(double *a, int size, double *x, double *y);

	// ��������� ������� ������� �� ������
	void Mult_MV_block_2x2(int nb, int *ig, int *jg, int *idi, int *ijg, 
		double *di_block, double *ggl_block, double *x, double *y, double *y_omp);

	void Mult_MV_block_2x2_transp(int nb, int *ig, int *jg, int *idi, int *ijg, 
		double *di_block, double *ggl_block, double *x, double *y, double *y_omp);


	///////////////// ��� ������-������������� ������������������ ////////////////////

	int Build_block_diag_preconditioner(int n, int *idi, double *di_block,
		double *df, int *idi_f,
		double *ggl_f, double *ggu_f);

	int Build_block_diag_preconditioner(int nb, int *idi, 
		double *di_block, double *df, double *ggl_f, double *ggu_f);

	int Build_complex_diag_preconditioner(int nb, int *idi, double *di_block, double *df);

	int solve_l_blockdiag(int n, double *df, int *idi_f, double *ggl_f, 
		double *f, double *x);

	int solve_u_blockdiag(int n, double *df, int *idi_f, double *ggl_f, 
		double *f, double *x);

	int solve_l_blockdiag(int n, double *df, double *ggl_f,double *f, double *x);
	int solve_u_blockdiag(int n, double *df, double *ggl_f, double *f, double *x);

	void mult_l_blockdiag(int nb, double *df, int *idi_f, double *ggl_f, double *x, double *y);
	void mult_l_blockdiag(int nb, double *df, double *ggl_f, double *x, double *y);

	void mult_u_blockdiag(int nb, double *df, int *idi_f, double *ggu_f, double *x, double *y);
	void mult_u_blockdiag(int nb, double *df, double *ggu_f, double *x, double *y);

	// ���������� ������������
	int LLT_Cmplx(int nb, int *ig, int *jg, int *idi, int *ijg, double *di, double *gg,
		double *d, double *sg);

	// ������� ���� � ���������������� ����������� ��������
	int SolveL_Cmplx(int nb, int *ig, int *jg, double *di, double *gg,
		double *f, double *x);

	// ������� ���� � ����������������� ����������� ��������
	int SolveU_Cmplx(int nb, int *ig, int *jg, double *di, double *gg,
		double *f, double *s, double *x);

	void Perm(int nb, double *x, double *y);
};

//------------------------------------------------------------------------