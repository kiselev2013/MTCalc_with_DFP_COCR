#pragma once
//--------------------------------------------
class Hex_Local_Matrix
{
public:
	int number_of_element;
	double x[8], y[8], z[8]; // ���������� ������ �������������	
	double J[3][3];          // ������� �����
	double J_1_T[3][3];		 // J^{-T}
	double det_J;            // ������� (������������ ������� �����)
	double det_J_abs;        // ������ ��������
	int type_of_hex; // 0..30 - ��������������, 31..61 - ������������
	double grad_all[8][3]; // ��������� �� ���. ������� � ����� ��������������
	// ��������� ������� 8�8 (�������, �� �������)
	double b[8][8]; // ������� �����
	double hx, hy, hz; // ��� ���������������
	bool JforParCalc; // ���������� �� ����� ������� � ������ ���������������
	void Calc_J(int n_of_point);
	void Calc_local_matrix_b_for_parallelepiped();
	void CalcMassMatrix();
	Hex_Local_Matrix(int i, long (*nver)[14], double (*xyz)[3]);
	~Hex_Local_Matrix();
};
