/**
 * GENERAL REMARKS
 *
 *  This code is freely available under the following conditions:
 *
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer.
 *
 *              MTCalc_with_DFP_COCR
 *  Header file for OutputResultant3d.cpp
 *
 *  Written by Prof. Marina G. Persova and Ph.D. Petr A. Domnikov 
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  mpersova@mail.ru
 *  Version 1.5 December 17, 2020
*/

#pragma once
#include "OutputArbitrary.h"
#include "Subdomain.h"
#include "AbstractFEM.h"

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

class OutputResultant3d
{
private:
	vector<Subdomain> sub;
	
	int n_pointres;
	double (*pointres)[3];

	AbstractFEM3D *TaskCalcMesh;

	int levelNeighbors;
	
	vector<long> ElemForPoint;
	vector< vector<long> > PointsForElem;
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
