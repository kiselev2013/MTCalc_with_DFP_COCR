#pragma once
class Divgrad_1d
{
public:
	double *coords_1d; 
	double *sin_1d;
	double *cos_1d;
	long n_1d; // ����� ����� � ���������� �����
	double bak; // ������ ������� ����
	double step0; // ��������� ��� (� ����������� �����) ��� ���������� �����
	double coef_razr; // ����������� ��������
	long n_layers_1d; // ����� ����
    double *layers_1d; // ������� ����
	double *sigma_1d;
	long alpha;

	int Solve_1d_Problem_for_3d_task();
	
	Divgrad_1d();
	~Divgrad_1d();
};
