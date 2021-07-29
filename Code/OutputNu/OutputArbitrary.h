#pragma once
#include <functional>
#include "T_Mapping.h"
#define OUTPUT
//------------------------------------------------------------------------
/*
	������� ������ �������� ��� �������� � ������������� ������������ ������
	���� ������������ �����������!!!!!!!!!!!!!!!!!!!

//  [4/2/2008 Domnikov]
	�������, ������ ������� ���������:

	- � ��� (�������) �������� �� ��� ���� ����� 2d-spline
		  - (������,
			��� �������, ��� ����������� �� x � �� y ������ ����� �������� �� 3-� ������ =>
			� ����� ������ ������������� ���������� �������� �� 1, 3(2 ��������), 9 ������)
		  - ��� ����� ����� �������� ����������� �� x,y ����� ������������� ����������� (2d),
		    � ��������� - ��������?????

	- ��� ��������� ���� � ��� (�������) - ������ ������ (� ����������)

	- � ����� (����):
		- ����� �������� ����������� �� �� �����������, �������� Bz ����� 3-������ ������
			(������������� �������� ������ ����� 1, 3 � 27 �����,
			�.�. Bz �� �������� ������� ������ �� z)
			- ����� � ���� ������� ������ �������� (�� ��� Bz ��������� ������������� �����������)
		- ���� �������� ����������� �� ����������� - �������� Bz ����� 2d-spline ��� ������ (1, 9 �����)

	- � ����� (����):
		- � ������, ���� �� ����������� �� ������� �������� - �� (Ex, Ey, Bz)
			������� ��� ���� ����� 2d-spline;
		- ���� �� ����������� ���� �������� � ����� ���������������� - ������� �� ����������
		- ���� �� ����������� ���� �������� � ����� �������������� - ������ ������ � �������������

	- � ����� (������� ���)
		- � ������, ���� �� ����������� �� ������� �������� - �� (Ex, Ey, Bz)
	������� ��� ���� ����� 2d-spline;
		- ���� �� ����������� ���� �������� � ����� ���������������� - ������� �� ����������
		- ���� �� ����������� ���� �������� � ����� �������������� - ������ ������ � �������������
	     
*/
//------------------------------------------------------------------------
// ������� � ������������ � ������� ��������� - ��� ����������
//------------------------------------------------------------------------
class PointRes2
{
public:
	long num; // ����� �����
	double point[3]; // ���������� �����

	PointRes2()	{}
	~PointRes2() {}
};
//------------------------------------------------------------------------
// ��������� ��� ���������� ���������
//------------------------------------------------------------------------
struct PointRes_less_x : public binary_function<PointRes2, PointRes2, bool>
{
	bool operator() (const PointRes2& a, const PointRes2& b) const {return a.point[0] < b.point[0];}
};
//------------------------------------------------------------------------
struct PointRes_less_y : public binary_function<PointRes2, PointRes2, bool>
{
	bool operator() (const PointRes2& a, const PointRes2& b) const {return a.point[1] < b.point[1];}
};
//------------------------------------------------------------------------
struct PointRes_less_z : public binary_function<PointRes2, PointRes2, bool>
{
	bool operator() (const PointRes2& a, const PointRes2& b) const {return a.point[2] < b.point[2];}
};
//------------------------------------------------------------------------
class Long_double_eq : public unary_function<long_double, bool>
{
	long_double a;
public:
	explicit Long_double_eq(const long_double& aa) : a(aa) {}
	bool operator() (const long_double& b) const {return fabs(b.d - a.d) < 1e-6;}
};
//------------------------------------------------------------------------
class Long_double_num_eq : public unary_function<long_double, bool>
{
	long_double a;
public:
	explicit Long_double_num_eq(const long_double& aa) : a(aa) {}
	bool operator() (const long_double& b) const {return a.i == b.i;}
};
//------------------------------------------------------------------------
struct Long_double_less : public binary_function<long_double, long_double, bool>
{
	bool operator() (const long_double& a, const long_double& b) const {return a.d < b.d;}
};
//------------------------------------------------------------------------
struct Long_double_num_less : public binary_function<long_double, long_double, bool>
{
	bool operator() (const long_double& a, const long_double& b) const {return a.i < b.i;}
};

//------------------------------------------------------------------------
// ������� ����� ��� ������
//------------------------------------------------------------------------
struct Output3dArbitrary
{
	int withSpline3d;
	int withSpline2d;
	int zeroPlane; // ������ �� ��������� z=0 ����� ����� ������ �� ������� ������ ��-���
	long kpar;
	long kuzlov;
	long n_pointres;
	double (*pointres)[3];
	double (*xyz)[3];

	long (*nver)[14];
	long (*nvtr)[8];
	long *type_of_hex; // 0..30 - ��������������, 31..61 - ������������

	vector<long> elemForPoint; // ��� ������� �������� ������ ��-�, ���� �� ��������
	vector< vector<long> > PointresForElem; // ��� ������� ��-�� ������ ������ ���������, �-��� � ���� ��������
	vector<long_double> PointresXsorted;
	vector<long_double> PointresYsorted;
	vector<long_double> PointresZsorted;
	
	Output3dArbitrary(int withSpline3d, int withSpline2d, int zeroPlane, long kuzlov, long kpar, long n_pointres,
		double (*pointres)[3], double (*xyz)[3]);
	~Output3dArbitrary();

	void PointresSort();
	int FindElemForReceivers();

	long GetGlobalVertex(long nElem, long nLocVertex);
	long GetTypeOfElement(long nElem);

	// ���������� ����� ���������, � ������� ��������  ����� ��� -1, ���� �� ��������
	int IsPointInsideHexahedron(double xm, double ym, double zm, double *x, double *y, double *z);

	int IsPointInsideBrick(double xm, double ym, double zm, double *x0, double *x1);

	// ����������������� ������� �� ������� (�� 3-������� �����) � ����� t �� [t_j2, t_j]
	double dA_dt(double t,
		double u_j, double u_j1, double u_j2, 
		double dt, double dt0, double dt1, 
		double t_j, double t_j1, double t_j2);
};

//------------------------------------------------------------------------
// ������ � ��������� ���������
//------------------------------------------------------------------------
struct OutputVect3d : public Output3dArbitrary
{
	long n_edges;
	long (*ed)[25];

	// ����������� ��� ��������� ������
	OutputVect3d(int withSpline3d, int withSpline2d, int zeroPlane,
		long kuzlov, long kpar, long n_edges, long n_pointres,
		double (*pointres)[3], double (*xyz)[3], long (*nver)[14], long (*ed)[25]);

	~OutputVect3d();

	int OutputFieldAtReceivers(double *v3, double *result, int derive);

	// ������ �������� �� ����� � ������ � ������������������� �� 3-������� �����
	int OutputDtAtReceivers(double *v3_j2, double *v3_j1, double *v3_j, double t,
		double t_j2, double t_j1, double t_j, int derive, double *result);
};

//------------------------------------------------------------------------
// ������ � ������� ���������
//------------------------------------------------------------------------
struct OutputNode3d : public Output3dArbitrary
{
	OutputNode3d(int withSpline3d, int withSpline2d, int zeroPlane,
		long kuzlov, long kpar, long n_pointres,
		double (*pointres)[3], double (*xyz)[3], long (*nvtr)[8], long *type_of_hex);

	~OutputNode3d();

	int OutputFieldAtReceiversHarm(double *v3, vector<vector<double>> &res, double w, int is_b);
	int OutputFieldAtReceivers(double *v3, vector<double> &result, int derive);
};
//------------------------------------------------------------------------

