#pragma once
class Global_slae_1d_harm_prof
{
public:

	long n;       // ����������� ����
	long n_elem;  // ����� ��������� � �����

	long *ig;    // ������ ������� ����� �����
	double *ggl; // �������� ������� ������������
	double *ggu; // �������� �������� ������������
	double *di;  // ���������
	double *pr;  // ������ ������ �����

	long ig_n_1; // ����������� ggl, ggu

	long *nvkat; // ������ ����������

	double *xyz; // �����

	long ay_hy;

	// ������������ ���������
	double *mu;
	double *sigma;
	double omega;

	Global_slae_1d_harm_prof(long n, long n_elem, long *ig, long *nvkat,
		double *mu, double *sigma, double omega, double *xyz, long ay_hy);
	~Global_slae_1d_harm_prof();

	void Assembling_for_1d_harm_problem();
	void Set_boundary_conditions_for_Ay();
	void Set_boundary_conditions_for_Hy();
};
