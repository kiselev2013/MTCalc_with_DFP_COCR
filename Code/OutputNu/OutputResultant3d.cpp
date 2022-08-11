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
 *  Functions for smoothing EM field
 *
 *  Written by Prof. Marina G. Persova and Ph.D. Petr A. Domnikov 
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  mpersova@mail.ru
 *  Version 1.5 December 17, 2020
*/
#include "stdafx.h"
#include "T_Brick.h"
#include "OutputResultant3d.h"
#include "Portret.h"
#include "los_rsf.h"
extern ofstream logfile;
//-----------------------------------------------------------
// Constructor
//-----------------------------------------------------------
OutputResultant3d::OutputResultant3d(AbstractFEM3D *pTaskCalcMesh,const Res3DValueType& r_type)
{
	ifstream inf;
	int i,n_neib;
	Maxiter=100;
	Nev=1e-40;
	pv::Point3D TmpPoint;
	TaskCalcMesh=pTaskCalcMesh;
	ValueType=r_type;
	n_pointres=TaskCalcMesh->GetNumberOfResPoints(ValueType);
	pointres = new double[n_pointres][3];
	for(i=0;i<n_pointres;i++){
		TmpPoint=TaskCalcMesh->GetResPoint(ValueType,i);
		pointres[i][0]=TmpPoint.x();
		pointres[i][1]=TmpPoint.y();
		pointres[i][2]=TmpPoint.z();
	}
	n_neib=5;
	inf.open("C:\\GI\\n_neib");
	if(inf)
	{
		inf>>n_neib;
		inf.close();
	}
	inf.clear();
	SetLevelNeighbors(n_neib);
}
//-----------------------------------------------------------
// Destructor
//-----------------------------------------------------------
OutputResultant3d::~OutputResultant3d()
{
	if(pointres) {delete [] pointres; pointres=NULL;}
}
//-----------------------------------------------------------
// Determine the number of girths to build the subdomain mesh
//-----------------------------------------------------------
void OutputResultant3d::SetLevelNeighbors(int val)
{
	levelNeighbors = val;
}
//-----------------------------------------------------------
// Âuild smoothing subdomains and define elements for output points
//-----------------------------------------------------------
int OutputResultant3d::Prepare()
{
	PointresSort();
	FindPointsForElems();
	InitSubdomains();
	return 0;
}
//-----------------------------------------------------------
// Âuild smoothing subdomains
//-----------------------------------------------------------
int OutputResultant3d::InitSubdomains()
{
	int i,j;
	vector<int> mat;
	int n_sub;

	int emat;
	int kpar=TaskCalcMesh->GetNumberOfElements();
	
	mat.clear();

	vRC.resize(n_pointres);

	for (i=0; i<kpar; i++)
	{
		int sz = (int)PointsForElem[i].size();

		if (sz > 0)
		{
			emat=TaskCalcMesh->GetElementMaterial(i);

			if (!binary_search(mat.begin(), mat.end(), emat)) 
			{
				mat.push_back(emat);
				sort(mat.begin(), mat.end());
			}
		
			double x[8],y[8],z[8];
			
			for(j=0;j<8;j++)
			{
				pv::Point3D p=TaskCalcMesh->GetNodeTrue(TaskCalcMesh->GetNodeNumberOnElement(i,j));
				x[j]=p.getx();
				y[j]=p.gety();
				z[j]=p.getz();
			}

			T_Brick L(x,y,z,TaskCalcMesh->GetTypeOfElement(i));
			for(j=0;j<sz;j++)
			{
				int pnt = PointsForElem[i][j];
				L.ScalarFieldOnParCff(pointres[pnt][0], pointres[pnt][1], pointres[pnt][2], vRC[pnt].cfa);
			}
		}
	}

	if(ValueType==vtWithoutDiscontinuity)
	{
		n_sub = 1;
		sub.resize(n_sub);
		sub[0].Init(-1, levelNeighbors, PointsForElem, TaskCalcMesh);
	}
	else
	{
		n_sub = (int)mat.size();
		sub.resize(n_sub);
		for (i=0; i<n_sub; i++)
		{
			sub[i].Init(mat[i], levelNeighbors, PointsForElem, TaskCalcMesh);
		}
	}

	return 0;
}
//-----------------------------------------------------------
// Output smooth values of EM fields
//-----------------------------------------------------------
int OutputResultant3d::Output(int itime)
{
	int i, j, k, l,n_sub;
	double var;

	n_sub=(int)sub.size();

	for (k=0; k<n_sub; k++)
	{
		Subdomain &mainSub=sub[k];

		for(j=0;j<mainSub.n_elem;j++)
		{
			mainSub.ValueInCenter[j]=
			TaskCalcMesh->GetValueInElementCenter(mainSub.renumElemFromNewToOld[j],
			ValueType);
		}
	
		mainSub.AsmGlobalVector();
	
		for (i=0; i<mainSub.n_nodes; i++){mainSub.x[i]=0.0;}

		mainSub.prds.solve_nrhs(1,mainSub.pr,mainSub.x);

		mainSub.CalcValuesAll_Node(mainSub.x);
		for (i=0; i<mainSub.n_elem; i++)
		{
			T_Brick L(i, mainSub.nver, mainSub.xyz);
			double ves[8];

			for (j=0; j<8; j++)
				ves[j] = mainSub.x[mainSub.nver[i][j]];

			int elem = mainSub.renumElemFromNewToOld[i];
			int sz = (int)PointsForElem[elem].size();

			for (j=0; j<sz; j++)
			{
				int pnt = PointsForElem[elem][j];
				var=0.0;
				for(l=0;l<8;l++){var+=ves[l]*vRC[pnt].cfa[l];}
				TaskCalcMesh->SaveResult(ValueType,var,pnt,itime);
			}
		}

	}

	return 0;
}
//-----------------------------------------------------------
// Sorting output points
//-----------------------------------------------------------
void OutputResultant3d::PointresSort()
{
	vector<PointRes2> p;
	long i;

	logfile << "xyz-Sorting of " << n_pointres << " receiver(s) start\n";

	p.resize(n_pointres);


	for (i=0; i<n_pointres; i++)
	{
		p[i].num = i;
		p[i].point[0] = pointres[i][0];
		p[i].point[1] = pointres[i][1];
		p[i].point[2] = pointres[i][2];
	}
	sort(p.begin(), p.end(), PointRes_less_x());

	PointresXsorted.resize(n_pointres);
	for (i=0; i<n_pointres; i++)
	{
		PointresXsorted[i].i = p[i].num;
		PointresXsorted[i].d = p[i].point[0];
	}
	logfile << "pointres x-sorting done"<< '\n';


	for (i=0; i<n_pointres; i++)
	{
		p[i].num = i;
		p[i].point[0] = pointres[i][0];
		p[i].point[1] = pointres[i][1];
		p[i].point[2] = pointres[i][2];
	}
	sort(p.begin(), p.end(), PointRes_less_y());
	PointresYsorted.resize(n_pointres);
	for (i=0; i<n_pointres; i++)
	{
		PointresYsorted[i].i = p[i].num;
		PointresYsorted[i].d = p[i].point[1];
	}
	logfile << "pointres y-sorting done"<< '\n';


	for (i=0; i<n_pointres; i++)
	{
		p[i].num = i;
		p[i].point[0] = pointres[i][0];
		p[i].point[1] = pointres[i][1];
		p[i].point[2] = pointres[i][2];
	}
	sort(p.begin(), p.end(), PointRes_less_z());
	PointresZsorted.resize(n_pointres);
	for (i=0; i<n_pointres; i++)
	{
		PointresZsorted[i].i = p[i].num;
		PointresZsorted[i].d = p[i].point[2];
	}
	logfile << "pointres z-sorting done"<< '\n';
}
//-----------------------------------------------------------
// Find elements for output points
//-----------------------------------------------------------
int OutputResultant3d::FindPointsForElems()
{
	long i, j, sz, t;
	long typeOfHex;

	vector<long_double>::iterator beg_x, beg_y, beg_z, end_x, end_y, end_z; 

	long_double coordMin[3], coordMax[3]; 

	vector<long_double> pntInBrick, pntInBrickTmp;
	vector<long_double> tmp_x, tmp_y, tmp_z;

	int kpar=TaskCalcMesh->GetNumberOfElements();

	logfile << "Start searching elements for receivers. ";
	logfile << "Elements=" << kpar << " receivers=" << n_pointres << endl;
	fflush(stdout);

	PointsForElem.resize(kpar);

	ElemForPoint.resize(n_pointres);
	for(i=0; i<n_pointres; i++)
		ElemForPoint[i] = -1;		

	double mdist;
	vector<int> edist;
	vector<double> vdist;

	edist.resize(n_pointres);
	vdist.resize(n_pointres);

	for(i=0;i<n_pointres;i++)
	{
		edist[i]=-1;
		vdist[i]=1e+30;
	}

	double recbrd[3][2];
	for(j=0;j<3;j++)
	{
		recbrd[j][0]=recbrd[j][1]=0.0;
	}
	if(n_pointres)
	{
		for(j=0;j<3;j++)
		{
			recbrd[j][0]=recbrd[j][1]=pointres[0][j];
		}
	}
	for(i=1;i<n_pointres;i++)
	{
		for(j=0;j<3;j++)
		{
			if(pointres[i][j]<recbrd[j][0])recbrd[j][0]=pointres[i][j];
			if(pointres[i][j]>recbrd[j][1])recbrd[j][1]=pointres[i][j];
		}
	}
	for(j=0;j<3;j++)
	{
		recbrd[j][0]-=1e-3;
		recbrd[j][1]+=1e-3;
	}

	for(i=0; i<kpar; i++)
	{
		typeOfHex = TaskCalcMesh->GetTypeOfElement(i);

		if(ValueType==vtWithDiscontinuity && TaskCalcMesh->GetElementMaterial(i)==0)
		{
			continue;
		}

		pv::Point3D TmpPoint;

		TmpPoint=TaskCalcMesh->GetNodeTrue(TaskCalcMesh->GetNodeNumberOnElement(i,0));
		coordMin[0].d=TmpPoint.x()-1e-6;
		coordMin[1].d=TmpPoint.y()-1e-6;
		coordMin[2].d=TmpPoint.z()-1e-6;
		TmpPoint=TaskCalcMesh->GetNodeTrue(TaskCalcMesh->GetNodeNumberOnElement(i,7));
		coordMax[0].d=TmpPoint.x()+1e-6;
		coordMax[1].d=TmpPoint.y()+1e-6;
		coordMax[2].d=TmpPoint.z()+1e-6;

		beg_x = lower_bound(PointresXsorted.begin(), PointresXsorted.end(), coordMin[0], Long_double_less());
		if(beg_x == PointresXsorted.end())
			continue;

		beg_y = lower_bound(PointresYsorted.begin(), PointresYsorted.end(), coordMin[1], Long_double_less());
		if(beg_y == PointresYsorted.end())
			continue;

		beg_z = lower_bound(PointresZsorted.begin(), PointresZsorted.end(), coordMin[2], Long_double_less());
		if(beg_z == PointresZsorted.end())
			continue;

		end_x = upper_bound(PointresXsorted.begin(), PointresXsorted.end(), coordMax[0], Long_double_less());
		if(end_x == PointresXsorted.begin())
			continue;

		end_y = upper_bound(PointresYsorted.begin(), PointresYsorted.end(), coordMax[1], Long_double_less());
		if(end_y == PointresYsorted.begin())
			continue;

		end_z = upper_bound(PointresZsorted.begin(), PointresZsorted.end(), coordMax[2], Long_double_less());
		if(end_z == PointresZsorted.begin())
			continue;

		if(distance(beg_x, end_x) <= 0)
			continue;

		if(distance(beg_y, end_y) <= 0)
			continue;

		if(distance(beg_z, end_z) <= 0)
			continue;

		tmp_x.clear();
		copy(beg_x, end_x, back_inserter(tmp_x));
		sort(tmp_x.begin(), tmp_x.end(), Long_double_num_less());

		tmp_y.clear();
		copy(beg_y, end_y, back_inserter(tmp_y));
		sort(tmp_y.begin(), tmp_y.end(), Long_double_num_less());

		tmp_z.clear();
		copy(beg_z, end_z, back_inserter(tmp_z));
		sort(tmp_z.begin(), tmp_z.end(), Long_double_num_less());

		set_intersection(tmp_x.begin(), tmp_x.end(), tmp_y.begin(), tmp_y.end(), back_inserter(pntInBrickTmp), Long_double_num_less()); 
		set_intersection(tmp_z.begin(), tmp_z.end(), pntInBrickTmp.begin(), pntInBrickTmp.end(), back_inserter(pntInBrick), Long_double_num_less()); 

		sz = (long)pntInBrick.size();
		for (j=0; j<sz; j++)
		{
			t = pntInBrick[j].i;

			if (ElemForPoint[t] == -1)
			{
				ElemForPoint[t] = i;
				PointsForElem[i].push_back(t);
			}
		}

		pntInBrick.clear();
		pntInBrickTmp.clear();

	}

	fflush(stdout);

	{
		char str[256],buf[256],resList[256];

		bool flag = false;

		for (t=0; t<n_pointres; t++)
		{
			if (ElemForPoint[t] == -1)
			{
				i=edist[t];
				if(i==-1)
				{
					sprintf(buf,"%d",t);
					if(strlen(resList)+strlen(buf)<256)
					{
						strcat(resList,buf);
					}
					flag = true;
				}
				else
				{
					ElemForPoint[t]=i;
					PointsForElem[i].push_back(t);
				}
			}
		}

		if (flag)
		{
			sprintf(str,"There are no elements for receiver(s) N");
			strcat(str,resList);
			logfile<<str<<endl;
			exit(1);
		}
	}

	return 0;
}
//-----------------------------------------------------------
// Finish pardiso in subdomains
//-----------------------------------------------------------
int OutputResultant3d::StopSolvers()
{
	int k, n_sub;
	int ret;
	n_sub=(int)sub.size();
	for(k=0;k<n_sub;k++)
	{
		Subdomain &mainSub=sub[k];
		mainSub.prds.stop_solver();
	}
	return 0;
}
