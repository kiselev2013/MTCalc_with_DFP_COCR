#pragma once
#include "OutputArbitrary.h"
#include "Subdomain.h"
#include "AbstractFEM.h"

struct _Plane_{
	pv::Vector N;
	double D;
	void set_D(pv::Point3D t){D=-(N.x()*t.x()+N.y()*t.y()+N.z()*t.z());}
};

struct ResCoef8
{
	double cfa[8];
	ResCoef8()
	{
		for(int i=0;i<8;i++)
		{
			cfa[i]=0.0;
		}
	}
};

//------------------------------------------------------------------------
class OutputResultant3d
{
private:
	vector<Subdomain> sub;					// ���������� �� ����������
	
	int n_pointres;
	double (*pointres)[3];

	AbstractFEM3D *TaskCalcMesh;			// ����� �� ������� ��������� ������

	int levelNeighbors;						// ������� �������� ��������� �������� � ����������
	
	vector<long> ElemForPoint;				// ��� ������� �������� ������ ����� ��-��, ���� �� ��������
	vector< vector<long> > PointsForElem;	// ��� ������� ��-�� ������ ������ ���������, �-��� � ���� ��������
	vector<long_double> PointresXsorted; 
	vector<long_double> PointresYsorted;
	vector<long_double> PointresZsorted;

	int InitSubdomains();
	void PointresSort();
	int FindPointsForElems();

public:
	OutputResultant3d(AbstractFEM3D *TaskCalcMesh,const Res3DValueType& r_type);
	~OutputResultant3d();

	int Prepare();
	int Output(int itime);

	Res3DValueType ValueType;

	void SetLevelNeighbors(int val);

	int Maxiter;
	double Nev;

	int StopSolvers();

	vector<ResCoef8> vRC;
};
//------------------------------------------------------------------------
