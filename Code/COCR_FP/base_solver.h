#pragma once
//------------------------------------------------------------------------
class Base_solver
{
// ����� Base_solver �������� ������������ ������������ ��� ���������� 
// ��������� (��������� ������������, ����� ������� � �.�.)
public:
	Base_solver();
	~Base_solver();

	inline double Scal(double *x, double *y, int n);     // ��������� ������������
	double Norm_Euclid(double *x, int n);         // ��������� ����� �������
	double Projection(double *vec, double *axis);  // �������� ������� �� ���
	double Relative_Error(double *analytic, double *numeric, int n); // ������������� �����������
	double Spline(double x, int n, double *xyz, double *values);

	//��������� ������� (� ������� �������) �� ������
	void Mult_Plot(double *a,double *pr,double *rez,int n); 

	// �������� �������
	int  Givens1(double& x, double& y, double& c, double& s);
	void Givens2(double& x, double& y, double c, double s);
	int  Givens(double *a, double *f, int n);

	// �������� ��� �� ������� �������
	int Undirect(double *a, double *b, double *x, int n);

	// ������� ���� � ���������� ��������, ������ ����������� �������
	// �������� ������ ���� ��������� ������������.
	// ��� ����������, ���� ������� ��������������� ��������
	int Solve_square_subdiag(double *a, double *b, double *x, int n);

	// ������ � ���� ������������� �������, � ������� �����, eps,
	// ����� �������� � ������� ������� ����
	int WriteKitChrono(char *fname, double residual, double eps, int iter, double time);
	int Write_kit(char *fname, double residual, double eps, int iter, int time);
	int Write_kit(char *fname, double residual, double eps, int iter, int time, int n);
	int Write_kit(char *fname, double residual, double eps, int iter, int time, int n, int size_jg);
	int Write_kit(char *fname, double residual, double eps, int iter, int time, double change_of_solution);


	// ���������� ����� �������� ��� ������� ���� � ���������������� (���)��������
	void SetPrcndItr(int n);

	// ������� ������������� ����� ������� � ��������� ��������
	double GetResidual();

	// ��������� ������������� ����� ������� � ��������� ��������
	void SetResidual(double residual);



	// ��� ���������������� ���������� (������� �������� � ������� �������, ����������� ����� (�������) - � std::complex<double>)

	// ��������� ������������ ��� ����������������� ��������
	std::complex<double> ScalCmplxTrue(double *x, double *y, int nb); 

	// ����������-����������� ��������� ������������
	std::complex<double> ScalCmplx(double *x, double *y, int nb); 

	// ��������� ������� �� ����������� �����
	void MultCmplxNumVect(std::complex<double> a, double *x, double *y, int nb);

	// ��������� ��������� ������ ������������ ������� �� ���������� ������� ������������ �������
	void MultCmplxVectVect(int nb, double *a, double *b, double *c);

	// ������� ��������� ������ ������������ ������� �� ���������� ������� ������������ �������
	void DivCmplxVectVect(int nb, double *a, double *b, double *c);

	// z = x + a*y
	void Cmplx_axpy(std::complex<double> a, double *x, double *y, double *z, int nb);

	double GetMinimiz(int n, double *ap, double *r);
};

//----------------------------------------------------
