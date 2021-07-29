#pragma once
#include "T_Mapping.h"
//---------------------------------------------------------------------------------
/*
����� Vec_Prep_Data ������ ��� �������� ����� � ������ ���������� � ������,
������� GeoPrep ����������� �� ���� � �����.

������� ���������� �������������� ��� ����� � ��� �� ��������� ��������������.
��� �� ��������� ���� �� ��������, �� ������� ��� ������ ����� � ����� ��������.

�� ����� 3dmeshregular ���� ����������� ������ x,y-����������. (z �� �����������)
*/
//---------------------------------------------------------------------------------
struct Tensor
{
	double val[3][3];

	Tensor()
	{
		Clear();
	}

	void Clear()
	{
		for(int i=0; i<3; i++)
		{
			for(int j=0; j<3; j++)
				val[i][j] = 0;
		}
	}

	Tensor& operator = (const double& d)
	{
		Clear();

		for(int j=0; j<3; j++)
			val[j][j] = d;

		return *this;
	}

	bool operator == (const double& d)
	{
		bool flag=true;
		int i, j;
		const double eps = 1e-6;

		for (i=0; i<3; i++)
		{
			for (j=0; j<3; j++)
			{
				if (i==j)
				{
					if(fabs(val[i][i] - d) > eps)
					{
						flag = false;
						break;
					}
				}
				else
				{
					if (val[i][j] != 0)
					{
						flag = false;
						break;
					}
				}
			}

			if (!flag)
				break;
		}

		return flag;
	}

	bool operator != (const double& d)
	{
		return !(*this == d);
	}
};
//------------------------------------------------------------------------

class Vec_Prep_Data
{
public:

	struct EnLine
	{
		long mtr;
		double es[3];
		double ec[3];
	};

	Vec_Prep_Data();
	~Vec_Prep_Data();

	int Read_prep_data(); // ������ ����� ��� ��� �� ���������
	int ReadPrepDataHarmLoop(char *pointres_fname); // ������ ����� ��� ����� � ������������� ����������
	int Read_mtz_1d(); // ��������� usin.dat, ucos.dat, alfa, nu (��� ��� ���)

	int Read_mesh_for_nonstat_problem(char *pointres_fname); // ��� ����� �� ���������
	int Read_infite0();  // ��������� ���-�� ����� � ������� �� infite.0

	// ������ 3dmeshregular, ��������� ���������� ����� � ����������=interval
	// (���� ����� ��������� ������ � ����� ������, ��� � 3dmeshregular)
	int Read_3dmeshregular(long interval); 

	int LoadReceivers(char *fname, int& _n_pointres, double (*(&_pointres))[3]);

	long n_pointres;
	double (*pointres)[3]; // ���������� ���������

	long maxiter; // ������������ ����� �������� ��� �������� (�� config)
	double eps;   // epsilon ��� �������� (�� config)
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
	double (*xyzt)[3];  // ���������� �����
	long n_layers_1d;  // ���-�� ��������� ���� (����������) ��� ��������� (�� sreda1d.ay)
	double *layers_1d; // ���� ��� ��������� (�� sreda1d.ay)
	double *sigma_1d;  // sigma0 ��� ��������� (�� sreda1d.ay)

	int LoadVectorE0ForLine(int n);
	vector< vector<EnLine> > EnForLine;

	int N_X,N_Y,N_Z;
	vector<int> regular;
	vector<double> Xcrd,Ycrd,Zcrd;

	// 3dmeshregular
	long n_mesh_regular_x;  // ����� ����� �� x � 3dmeshregular
	long n_mesh_regular_y;  // ����� ����� �� y � 3dmeshregular
	double *mesh_regular_x; // x-���������� �� 3dmeshregular
	double *mesh_regular_y; // y-���������� �� 3dmeshregular

	// ������ ��� �������������� ������
	long ntime;     // ����� ��������� ����
	double *time;   // ��������� ����

	// ������ ��� ��� 
 	//long norvect; // ����� �������� ��� GMRES
 	double nu;    // �������
 	long alfa;    // ����������� ����: �� x - alpha=1; �� y - alpha=0; J=(alpha*Jx,(1-alpha)*Jy, 0)
 	double *usin; // ������� ���������� ������ (sin-����������)
 	double *ucos; // ������� ���������� ������ (sin-����������)
 	double *z_1d; // ���������� �����
 	long n_1d;    //  ����� ����� � ���������� �����

	T_Mapping_Vec *tmap;
	long tasktype;

	int npr, nfreq;

	int fdirect;
};
//-----------------------------------------------------------