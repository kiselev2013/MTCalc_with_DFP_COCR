#pragma once
class Local_matrix_1d
{
public:

	double b[2][2]; // ������� ��������
	double c[2][2]; // ������� ����� 
	double f_re[2], f_im[2]; // �������� ������ ����� � �����
	double g_re[2], g_im[2]; // ������� ������ ����� ��� Re � Im ���������
	double a[4][4]; // ������� ��� ������������� ������
	double g[4]; // ������ ������ ����� ��� ������������� ������

	// ������������ ���������
	double mu;
	double sigma;
	double omega;

	double h; // ������ ��������� ��������

	long alpha; // ����������� ����


	Local_matrix_1d(double h, double mu, double sigma, double omega,
		double *f_re, double *f_im, long alpha);
	~Local_matrix_1d();

	void Compute_local_matrix_and_vector_harm();

	void Compute_local_matrix_harm();
	void Compute_local_vector_harm();
	void Compute_local_matrix_b();
	void Compute_local_matrix_c();
	void Compute_local_vector(double *f, double *g);
};
