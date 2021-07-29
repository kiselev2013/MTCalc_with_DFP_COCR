#pragma once

class Vec_Prep_Data
{
public:

	Vec_Prep_Data();
	~Vec_Prep_Data();

	int Read_prep_data(); // ������ ����� ��� ��� �� ���������
	int Read_mtz_1d(); // ��������� usin.dat, ucos.dat, alfa, nu (��� ��� ���)

	int Read_3dmeshregular(long interval); 

	int LoadReceivers(char *fname, int& _n_pointres, double (*(&_pointres))[3]);

	long n_materials; // ����� ��������� ����������
	double *mu3d;     //(mu ��� ��������� ������ � ��������� (������ ������� �� ����� mu3d))
	double *mu0;      //(mu ��� ��������� (������ ������� �� ����� mu3d))
	int n_pointresB;  // ����� ��������� ��� B
	int n_pointresE;  // ����� ��������� ��� E
	double (*pointresB)[3]; // ���������� ���������
	double (*pointresE)[3]; // ���������� ���������
	double *sigma3d;       //(sigma ��� ��������� ������ � ��������� (������ ������� �� sig3d))
	double *sigma0;        //(sigma ��� ��������� (������ ������� �� sig3d))
	double *dpr3d;       //(dpr ��� ��������� ������ � ��������� (������ ������� �� dpr3d))
	double *dpr0;        //(dpr ��� ��������� (������ ������� �� dpr3d))
	long kuzlov;      // ����� ����� (����, � ��� ����� � ������������)
	long kpar;        // ����� ��������� � �����
	long kt1;         // ����� ����� � ������� ��������
	long *l13d;       // ������ ����� � ������� ��������
	long (*nver)[14]; // ������ ����� �������� ��������� (13-������� �������������)
	long *nvkat;       // ������ ���������� �������� ��-���
	double (*xyz)[3];  // ���������� �����
	long n_layers_1d;  // ���-�� ��������� ���� (����������) ��� ��������� (�� sreda1d.ay)
	double *layers_1d; // ���� ��� ��������� (�� sreda1d.ay)
	double *sigma_1d;  // sigma0 ��� ��������� (�� sreda1d.ay)

	int N_X,N_Y,N_Z;
	vector<int> regular;
	vector<double> Xcrd,Ycrd,Zcrd;

	long n_mesh_regular_x;  // ����� ����� �� x � 3dmeshregular
	long n_mesh_regular_y;  // ����� ����� �� y � 3dmeshregular
	double *mesh_regular_x; // x-���������� �� 3dmeshregular
	double *mesh_regular_y; // y-���������� �� 3dmeshregular

 	double nu;    // �������
 	long alfa;    // ����������� ����: �� x - alpha=1; �� y - alpha=0; J=(alpha*Jx,(1-alpha)*Jy, 0)
 	double *usin; // ������� ���������� ������ (sin-����������)
 	double *ucos; // ������� ���������� ������ (sin-����������)
 	double *z_1d; // ���������� �����
 	long n_1d;    //  ����� ����� � ���������� �����

	int npr, nfreq;
};
//-----------------------------------------------------------