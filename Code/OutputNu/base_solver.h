#pragma once

class Base_solver
{
/*
����� Base_solver �������� ������������ ������������ ��� ���������� 
��������� (��������� ������������, ����� ������� � �.�.)
*/

// �������� � ��������� ���������� � 1

protected:
	double *y_omp;


public:
	Base_solver();
	~Base_solver();

	int n_threads; // ����� ������� (����������� ��� OpenMP)

	double Scal(double *x, double *y, long n);     // ��������� ������������
	double Norm_Euclid(double *x, long n);         // ��������� ����� �������
	double Projection(double *vec, double *axis);  // �������� ������� �� ���
	double Relative_Error(double *analytic, double *numeric, long n); // ������������� �����������
	double Spline(double x, long n, double *xyz, double *values);

	//��������� ������� (� ������� �������) �� ������
	void Mult_Plot(double *a,double *pr,double *rez,long n); 

	// �������� �������
	int  Givens1(double& x, double& y, double& c, double& s);
	void Givens2(double& x, double& y, double c, double s);
	int  Givens(double *a, double *f, long n);

	// �������� ��� �� ������� �������
	int Undirect(double *a, double *b, double *x, long n);

	// ������� ���� � ���������� ��������, ������ ����������� �������
	// �������� ������ ���� ��������� ������������.
	// ��� ����������, ���� ������� ��������������� ��������
	int Solve_square_subdiag(double *a, double *b, double *x, long n);

	// ������ � ���� ������������� �������, � ������� �����, eps,
	// ����� �������� � ������� ������� ����
	int Write_kit(char *fname, double residual, double eps, long iter, __time64_t time);
	int Write_kit(char *fname, double residual, double eps, long iter, __time64_t time, long n);
	int Write_kit(char *fname, double residual, double eps, long iter, __time64_t time, long n, long size_jg);
	int Write_kit(char *fname, double residual, double eps, long iter, __time64_t time, double change_of_solution);


	// ��� ���������������� ���������� (������� �������� � ������� �������, ����������� ����� (�������) - � std::complex<double>)

	// ��������� ������������ ��� ����������������� ��������
	std::complex<double> ScalCmplx(double *x, double *y, long nb); 

	// ��������� ������� �� ����������� �����
	void MultCmplxNumVect(std::complex<double> a, double *x, double *y, long nb);

	// ��������� ��������� ������ ������������ ������� �� ���������� ������� ������������ �������
	void MultCmplxVectVect(long nb, double *a, double *b, double *c);

	// ������� ��������� ������ ������������ ������� �� ���������� ������� ������������ �������
	void DivCmplxVectVect(long nb, double *a, double *b, double *c);

protected:
	long nMRSrestart; // ����� ��������, ����� ������� ��������������� MRS
	// (����� �� ������������� ������ ����������)
public:
	void Set_nMRSrestart(long nMRSrestart) {this->nMRSrestart = nMRSrestart;}

};
//----------------------------------------------------
