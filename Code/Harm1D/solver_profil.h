#pragma once
class Solver_profil
{
public:
	long n;
	long *ig;
	double *di;
	double *ggl;
	double *ggu;
	double *pr;

	double *x;

    double zero; // ��� ��������� � ����

	double *h1, *h2; //��������������� �������
	bool mem_h; // ���������� �� ������ ��� ��������������� ������

	Solver_profil();
	Solver_profil(long n, long *ig, double *di, double *ggl ,double *ggu, double *pr, double *x);
	~Solver_profil();


	int LU_profil();
	int Solve_SLAE_using_LU();
	int Solve_SLAE_using_LU(long n, long *ig, double *di, double *ggl ,double *ggu, double *pr, double *x);

	int LU_profil(long n, long *ig, double *di, double *ggl ,double *ggu, double *sl, double *su, double *d);
	int LLT_profil(long n, long *ig, double *di, double *gg, double *sl, double *d);

	double Scal(double *a, double *b, long n);

	// ������� ���� � ���������������� �������� (1 �� ���������)
	void Solve_L_1(long n, long *ig, double *ggl, double *pr, double *y);

	// ������� ���� � ����������������� ��������
	int Solve_U(long n, long *ig, double *di, double *ggu, double *y, double *x, double *h);

	// ������� ���� � ���������������� ��������
	int Solve_L(long n, long *ig, double *di, double *ggl, double *y, double *x);

	int Solve_SLAE_using_LLT();

	int Solve_SLAE_using_LLT(long n, long *ig, double *di, double *ggl, double *pr, double *x);


};
