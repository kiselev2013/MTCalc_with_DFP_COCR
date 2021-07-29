#pragma once
#include <functional>
#include "T_Mapping.h"
#define OUTPUT
//------------------------------------------------------------------------
/*
	Вариант выдачи напрямую без сплайнов и согласованных результантов должен
	быть предусмотрен обязательно!!!!!!!!!!!!!!!!!!!

//  [4/2/2008 Domnikov]
	Вообщем, решено сделать следующее:

	- в МТЗ (узловом) оставить всё как есть через 2d-spline
		  - (однако,
			мне кажется, что производные по x и по y точнее всего выдавать по 3-м точкам =>
			в новой выдаче предусмотреть построение сплайнов по 1, 3(2 варианта), 9 точкам)
		  - или лучше всего выдавать производные по x,y через согласованный интерполянт (2d),
		    а потенциал - напрямую?????

	- для выделения поля в МТЗ (узловом) - прямая выдача (с тетраэдров)

	- в петле (ВМКЭ):
		- когда приёмники расположены не на поверхности, выдавать Bz через 3-мерный сплайн
			(предусматреть варианты выдачи через 1, 3 и 27 точек,
			т.к. Bz на элементе зависит только от z)
			- также я хочу сделать выдачу напрямую (но для Bz построить согласованный интерполянт)
		- если приёмники расположены на поверхности - выдавать Bz через 2d-spline как обычно (1, 9 точек)

	- в линии (ВМКЭ):
		- в случае, если на поверхность не выходят аномалии - всё (Ex, Ey, Bz)
			выдаётся как есть через 2d-spline;
		- если на поверхности есть аномалии в форме параллелепипедов - разбить на подобласти
		- если на поверхности есть аномалии в форме шестигранников - прямая выдача с треугольников

	- в линии (узловой МКЭ)
		- в случае, если на поверхность не выходят аномалии - всё (Ex, Ey, Bz)
	выдаётся как есть через 2d-spline;
		- если на поверхности есть аномалии в форме параллелепипедов - разбить на подобласти
		- если на поверхности есть аномалии в форме шестигранников - прямая выдача с треугольников
	     
*/
//------------------------------------------------------------------------
// приёмник с координатами и номером приемника - для сортировки
//------------------------------------------------------------------------
class PointRes2
{
public:
	long num; // номер точки
	double point[3]; // координаты точки

	PointRes2()	{}
	~PointRes2() {}
};
//------------------------------------------------------------------------
// предикаты для сортировки приёмников
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
// Базовый класс для выдачи
//------------------------------------------------------------------------
struct Output3dArbitrary
{
	int withSpline3d;
	int withSpline2d;
	int zeroPlane; // выдача на плоскости z=0 через точки Гаусса на верхних гранях эл-тов
	long kpar;
	long kuzlov;
	long n_pointres;
	double (*pointres)[3];
	double (*xyz)[3];

	long (*nver)[14];
	long (*nvtr)[8];
	long *type_of_hex; // 0..30 - параллелепипед, 31..61 - шестигранник

	vector<long> elemForPoint; // для каждого приёмника хранит эл-т, куда он попадает
	vector< vector<long> > PointresForElem; // для каждого эл-та хранит номера приёмников, к-рые в него попадают
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

	// возвращает номер тетраэдра, в который попадает  точка или -1, если не попадает
	int IsPointInsideHexahedron(double xm, double ym, double zm, double *x, double *y, double *z);

	int IsPointInsideBrick(double xm, double ym, double zm, double *x0, double *x1);

	// дифференцирование решения по времени (по 3-слойной схеме) в точке t из [t_j2, t_j]
	double dA_dt(double t,
		double u_j, double u_j1, double u_j2, 
		double dt, double dt0, double dt1, 
		double t_j, double t_j1, double t_j2);
};

//------------------------------------------------------------------------
// выдача с векторных элементов
//------------------------------------------------------------------------
struct OutputVect3d : public Output3dArbitrary
{
	long n_edges;
	long (*ed)[25];

	// конструктор для векторной задачи
	OutputVect3d(int withSpline3d, int withSpline2d, int zeroPlane,
		long kuzlov, long kpar, long n_edges, long n_pointres,
		double (*pointres)[3], double (*xyz)[3], long (*nver)[14], long (*ed)[25]);

	~OutputVect3d();

	int OutputFieldAtReceivers(double *v3, double *result, int derive);

	// выдать значения из сетки в точках и продифференцировать по 3-слойной схеме
	int OutputDtAtReceivers(double *v3_j2, double *v3_j1, double *v3_j, double t,
		double t_j2, double t_j1, double t_j, int derive, double *result);
};

//------------------------------------------------------------------------
// выдача с узловых элементов
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

