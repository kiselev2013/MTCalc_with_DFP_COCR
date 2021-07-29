#pragma once
#include "Vec_Prep_Data.h"

class T_Brick
{
public:
	double hx, hy, hz; // ������� ��-��
	double xk, xk1, yk, yk1, zk, zk1; // ������� ���������������
    long num; // ����� ��������� �������� � �����

	double b[12][12]; // ��������� ������� ��������
	double c[12][12]; // ��������� ������� �����

	double f_re[12];  // �������� ������ ����� � ��������� ����

	double a[12][12]; // ��������� ������� ��-��  
	double g[12];     // ��������� ������ ��-��
	double g_harm[24];
	double g_re[12], g_im[12];
	double g_re_b[12], g_im_b[12];

	// ������������ ��������� 
	double mu;
	double mu0;
	double sigma;
	double sigma0;
	double dpr;
	double dpr0;
	long n_mat; // ����� ��������� ��������

	long (*nver)[14]; // ������ ����� �������� ��������� (13-������� �������������)
	long (*ed)[25];   // ��-�� ������������� ������ ������ + ������������ ���� + ��� ��-��
	long *nvkat;      // ������ ���������� �������� ��-���
	long (*edges)[2]; // ����, �������� 2-�� ���������
	double (*xyz)[3]; // ���������� �����

	T_Brick(double *x_coords, double *y_coords, double *z_coords);

	T_Brick(long num, long (*nver)[14], long (*ed)[25], long (*edges)[2], double (*xyz)[3], long *nvkat,
		double *sigma3d, double *sigma0, double *mu3d, double omega,
		long alpha, long n_1d, double *z_1d, double *sin_1d, double *cos_1d);

	T_Brick(long num, long (*nver)[14], double (*xyz)[3]);

	~T_Brick();

	void Compute_Local_Matrix_And_Vector(const long what_compute); 
	void Compute_Local_Matrix_B(); 
	void Compute_Local_Matrix_C(); 

	double omega; // ����������� ������� ��� ������������� ������
	double asin0[12], acos0[12]; // ���������� ���� � ��������� ���� 

	long alpha;
	long n_1d;
	double *z_1d;
	double *sin_1d;
	double *cos_1d;

	void Calc_local_vector_for_MT();

	void Calc_asin_acos_at_middle_of_edges();

	void Calc_block_local_matrix_and_vector();

	double x[8], y[8], z[8]; // ���������� ������ �������������	
	double J[3][3];          // ������� �����
	double J_1[3][3];		 // J^{-1}
	double J_1_T[3][3];		 // J^{-T}
	double det_J;            // ������� (������������ ������� �����)
	double det_J_abs;        // ������ ��������
	double phi_all[12][3];    // �������� ������� � ����� ��������������
	double rot_all[12][3];    // ������ �� ���. ������� � ����� ��������������

	void Mapping(double *in, double *out);

	void Calc_J(int n_of_point); // ��������� ������� ����� � ����� ������ 
	void Calc_J(double x, double y, double z); // ��������� ������� ����� � ������������ ����� 
	void Calc_J_on_face(int n_of_point); // ��������� ������� � ������, ������������� � �����
	void Calc_J_in_parallelepiped(); // ��������� ������� ����� � ������ ���������������

	// ��������� �������� ���������� ���� ������ ������������� �������������
	void Calc_value_inside_hex(double *ves, double *in, double *out); 

	// �������� i-� �������� ������� �� ��������������� 
	void Basis_func_on_vec_par(long i, double *in, double *out);

	// �������� i-� �������� ������� �� ��������������� � ��������������� ������������� 2/hx, 2/hy, 2/hz
	void Basis_func_on_vec_par(long i, double ves, double *in, double *out);

	// �������� i-� �������� ������� �� ��������� ��������������� 
	void Basis_func_on_reference_vec_par(long i, double *in, double *out);	

	// �������� i-� �������� ������� �� ������������� 
	void Basis_func_on_vec_hex(long i, double ves, double *in, double *out);

	// ������ ����� � ������������ ����� ������ �������������
	void Calc_rotor_inside_hex(double *ves, double *in, double *out); 

	// �������� ������ i-� �������� ������� �� ��������� ��������������� 
	void Rot_of_basis_func_on_reference_vec_par(long i, double *in, double *out);

	// �������� ������ x-���������� ������ i-� �������� ������� ������ ���������������
	void Rotx_of_basis_func_on_reference_vec_par(long i, double x, double *out);

	// �������� ������ y-���������� ������ i-� �������� ������� ������ ���������������
	void Roty_of_basis_func_on_reference_vec_par(long i, double y, double *out);

	// �������� ������ z-���������� ������ i-� �������� ������� ������ ���������������
	void Rotz_of_basis_func_on_reference_vec_par(long i, double z, double *out);

	// �������� ������ i-� �������� ������� ������ ���������������
	void Rot_of_basis_func_on_vec_par(long i, double ves, double *in, double *out);

	// ������ z-���������� ������ � ������ ������ � ������� ����� �������������
	// ����� ��� ��� ������� (��� 3-������� ����� �� �������)
	void Get_rotz_on_face(double *ves1, double *ves2, double *ves3,
		double *out1, double *out2, double *out3);

	// ��������� ���������� �������� ������� �� [-1, 1]^3
	double l0(double x);
	double l1(double x);

	// ��������� ������� �������� ������� �� [-1, 1]^3
	double Phi_node(long i, double x, double y, double z);

	// ����������� �� ������� ��������� �������� �-���
	double dPhi_node(long i, long j, double x, double y, double z);

	// ������ �������� ���������� ���� �� ���������������
	// (���������� ����� - ����������)
	void VectorFieldOnPar(double x, double y, double z, double *ves,
		double *x_out, double *y_out, double *z_out);
	double ScalarFieldOnPar(double x, double y, double z, double *ves);
	double DxOfScalarFieldOnPar(double x, double y, double z, double *ves);
	double DyOfScalarFieldOnPar(double x, double y, double z, double *ves);
	double DzOfScalarFieldOnPar(double x, double y, double z, double *ves);

	// ������ x-���������� ���� � ��������������� ����� ��� ��� �������� �����
	void VectorFieldXOnPar3(double y, double z, double *ves_j2, double *ves_j1, double *ves_j,
		double *out_j2, double *out_j1, double *out_j);

	// ������ y-���������� ���� � ��������������� ����� ��� ��� �������� �����
	void VectorFieldYOnPar3(double x, double z, double *ves_j2, double *ves_j1, double *ves_j,
		double *out_j2, double *out_j1, double *out_j);

	void RotXOnPar(double x, double *ves, double *out, bool loc_c=false);
	void RotYOnPar(double y, double *ves, double *out, bool loc_c=false);
	void RotZOnPar(double z, double *ves, double *out);
	void RotZOnPar3(double z,
		double *ves1, double *ves2, double *ves3,
		double *out1, double *out2, double *out3);

	// ������ ����������
	void Transformation_of_variables(double *in, double *out);
	void Transformation_of_variables(double *x, double *y, double *z);
	double Xi(double x);
	double Eta(double y);
	double Zeta(double z);

	void GetAonFace(double *ves1, double *ves2, double *ves3,
		double *out1, double *out2, double *out3);

	void Set_dpr(double dpr);
	void Set_dpr0(double dpr0);
	void Set_mu0(double mu0);
	Vec_Prep_Data *d;

	double asin0n[8][3], acos0n[8][3]; // ���������� ���� � �����
	double asin0c[3], acos0c[3]; // ���������� ���� � ������ ��������

	void Calc_asin_acos_at_nodes();

	void GetVectorFieldNodes(double *ves, double *ax, double *ay, double *az);

	int n_edges;
};
//------------------------------------------------------------------------------
const double MIDDLE_OF_LOCAL_EDGE[12][3] = {
	 0.0, -1.0, -1.0,
	 0.0,  1.0, -1.0,
	 0.0, -1.0,  1.0,
	 0.0,  1.0,  1.0,

	-1.0,  0.0, -1.0,
	-1.0,  0.0,  1.0,
	 1.0,  0.0, -1.0,
	 1.0,  0.0,  1.0,

	-1.0, -1.0,  0.0,
	 1.0, -1.0,  0.0,
	-1.0,  1.0,  0.0,
	 1.0,  1.0,  0.0
};
//-----------------------------------------------------------------------
const double LOCAL_COORDS_OF_NODES[8][3] = 
{
	-1.0, -1.0, -1.0,
	 1.0, -1.0, -1.0,
	-1.0,  1.0, -1.0,
	 1.0,  1.0, -1.0,

	-1.0, -1.0,  1.0,
	 1.0, -1.0,  1.0,
	-1.0,  1.0,  1.0,
	 1.0,  1.0,  1.0
};
//-----------------------------------------------------------------------
const double TANGENT_VECTORS_ON_REFERENCE_CUBE[12][3] = {
	1.0, 0.0, 0.0,
	1.0, 0.0, 0.0,
	1.0, 0.0, 0.0,
	1.0, 0.0, 0.0,

	0.0, 1.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 1.0, 0.0,

	0.0, 0.0, 1.0,
	0.0, 0.0, 1.0,
	0.0, 0.0, 1.0,
	0.0, 0.0, 1.0
};
//-----------------------------------------------------------------------
const long REG_EDGES[12][2]={ 	// ����� �������������� ���� ���� �� ��������
	0,1, 2,3, 4,5, 6,7, 0,2, 4,6, 1,3, 5,7, 0,4, 1,5, 2,6, 3,7 };
//-----------------------------------------------------------------------