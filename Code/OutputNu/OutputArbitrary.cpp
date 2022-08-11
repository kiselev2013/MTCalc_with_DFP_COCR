/**                                                                                       
 * GENERAL REMARKS                                                                        
 *                                                                                        
 *  This code is freely available under the following conditions:                         
 *                                                                                        
 *  1) The code is to be used only for non-commercial purposes.                           
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer. 
 *                                                                                        
 *  			MTCalc_with_DFP_COCR                                              
 *  This file contains some basic routines for output of values of the desired function   
 *  at an arbitrary point in the computational domain (3D spline)                         
 *                                                                                        
 *  Written by Ph.D. Petr A. Domnikov                                                     
 *  Novosibirsk State Technical University,                                               
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                     
 *  p_domnikov@mail.ru                                                                    
 *  Version 1.3 January 10, 2021                                                          
*/                                                                                        

#include "stdafx.h"
#include "T_Brick.h"
#include "OutputArbitrary.h"
extern ofstream logfile;
//------------------------------------------------------------------------
Output3dArbitrary::Output3dArbitrary(int withSpline3d, int withSpline2d, int zeroPlane,
									 long kuzlov, long kpar, long n_pointres,
									 double (*pointres)[3], double (*xyz)[3])
{
	this->withSpline3d = withSpline3d;	
	this->withSpline2d = withSpline2d;	
	this->zeroPlane = zeroPlane;
	this->kuzlov = kuzlov;
	this->kpar = kpar;
	this->n_pointres = n_pointres;
	this->pointres = pointres;
	this->xyz = xyz;

	this->nver = NULL;
	this->nvtr = NULL;
	this->type_of_hex = NULL;

	PointresSort();
}
//------------------------------------------------------------------------
Output3dArbitrary::~Output3dArbitrary()
{
}
//---------------------------------------------------------------------------------------
// differentiation of the solution with respect to time (according to the 3-layer scheme) at the point t from [t_j2, t_j]
//---------------------------------------------------------------------------------------
double Output3dArbitrary::dA_dt(double t, double u_j, double u_j1, double u_j2,
								double dt, double dt0, double dt1, double t_j, double t_j1, double t_j2)
{
	double du_dt;

	if(t < t_j2) t = t_j2;
	if(t > t_j)  t = t_j;

	du_dt = u_j2*(2.0*t - t_j  - t_j1)/(dt1*dt) - 
		u_j1*(2.0*t - t_j  - t_j2)/(dt1*dt0) + 
		u_j*(2.0*t - t_j1 - t_j2)/(dt*dt0);

	return du_dt;
}
//------------------------------------------------------------------------
// constructor for vector problem
//------------------------------------------------------------------------
OutputVect3d::OutputVect3d(int withSpline3d, int withSpline2d, int zeroPlane, long kuzlov, long kpar, long n_edges,
						   long n_pointres, double (*pointres)[3], double (*xyz)[3], long (*nver)[14],
						   long (*ed)[25]) : Output3dArbitrary(withSpline3d, withSpline2d, zeroPlane,
						   kuzlov, kpar, n_pointres, pointres, xyz)
{
	this->n_edges = n_edges;
	this->nver = nver;
	this->ed = ed;

	FindElemForReceivers();
}
//------------------------------------------------------------------------
OutputVect3d::~OutputVect3d()
{
}
//------------------------------------------------------------------------
// constructor for node problem
//------------------------------------------------------------------------
OutputNode3d::OutputNode3d(int withSpline3d, int withSpline2d, int zeroPlane,
						   long kuzlov, long kpar, long n_pointres, double (*pointres)[3],
						   double (*xyz)[3], long (*nvtr)[8], long *type_of_hex) : Output3dArbitrary(withSpline3d,
						   withSpline2d, zeroPlane, kuzlov, kpar, n_pointres, pointres, xyz)
{
	this->nvtr = nvtr;
	this->type_of_hex = type_of_hex;
	
	FindElemForReceivers();
}
//------------------------------------------------------------------------
OutputNode3d::~OutputNode3d()
{
}
//------------------------------------------------------------------------
// output values from the grid in points and differentiate according to the 3-layer scheme
// for output of E=-dA/at and Eds=-dBz/dt
//------------------------------------------------------------------------
int OutputVect3d::OutputDtAtReceivers(double *v3_j2, double *v3_j1, double *v3_j, double t,
						double t_j2, double t_j1, double t_j, int derive, double *result)
{
	const double dt =  t_j  - t_j2;
	const double dt0 = t_j  - t_j1;
	const double dt1 = t_j1 - t_j2;

	long i, j;
	long sz;
	long nElem;
	long typeOfElem;
	double x[8], y[8], z[8]; // coordinates of the current hexahedron, from which we output the result
	double tx, ty, tz; // current receiver coordinates
	double ves_j[12], ves_j1[12], ves_j2[12]; // weights for the edges of the current hexahedron on three time layers
	double val_j2, val_j1, val_j;

	if(withSpline3d)
	{
		if (withSpline2d)
		{
		} 
		else
		{
		}
	}
	else
	{
		if (withSpline2d)
		{
		} 
		else
		{
			for (nElem=0; nElem<kpar; nElem++)
			{
				sz = (long)PointresForElem[nElem].size(); // the number of receivers that are in this element
				if(sz == 0)
					continue;

				typeOfElem = GetTypeOfElement(nElem);

				if (typeOfElem <= 30) // parallelepiped
				{
					// initialization of an object of class T_Brick
					for (i=0; i<8; i++)
					{
						j = GetGlobalVertex(nElem, i);
						x[i] = xyz[j][0];
						y[i] = xyz[j][1];
						z[i] = xyz[j][2];
					}
					T_Brick L(x, y, z, typeOfElem);

					// determine the weights for the basis functions
					for (i=0; i<12; i++)
					{
						ves_j2[i] = v3_j2[ed[nElem][i]];
						ves_j1[i] = v3_j1[ed[nElem][i]];
						ves_j[i]  = v3_j[ed[nElem][i]];
					}

					// for each point that falls into the element, it outputs the required value from the parallelepiped 
					for (i=0; i<sz; i++)
					{
						j = PointresForElem[nElem][i];
						tx = pointres[j][0];
						ty = pointres[j][1];
						tz = pointres[j][2];

						// output values on three time layers
						switch(derive)
						{
						case 0: // Ex = -dAx/dt
							L.VectorFieldXOnPar3(ty, tz, ves_j2, ves_j1, ves_j, &val_j2, &val_j1, &val_j);
							break;
						case 1: // Ey = -dAy/dt
							L.VectorFieldYOnPar3(tx, tz, ves_j2, ves_j1, ves_j, &val_j2, &val_j1, &val_j);
							break;
						case -3: // Eds = -dBz/dt
							L.RotZOnPar3(tz, ves_j2, ves_j1, ves_j, &val_j2, &val_j1, &val_j);
							break;
						}
						
						// differentiate with respect to time
						result[j] = -dA_dt(t, val_j, val_j1, val_j2, dt, dt0, dt1, t_j, t_j1, t_j2);
					}
				} 
				else // hexahedron
				{
				}
			}				
		}
	}

	return 0;
}

//------------------------------------------------------------------------
int OutputVect3d::OutputFieldAtReceivers(double *v3, double *result, int derive)
{
	long i;
	long t;
	long sz;
	long nElem;
	long typeOfElem;
	double x[8], y[8], z[8];
	double tx, ty, tz;
	double ves[12];
	double tmp1, tmp2;

	if(withSpline3d)
	{
		if (withSpline2d)
		{
		} 
		else
		{
		}
	}
	else
	{
		if (withSpline2d)
		{
		} 
		else
		{
			for (nElem=0; nElem<kpar; nElem++)
			{
				sz = (long)PointresForElem[nElem].size(); // the number of receivers that are in the element
				if(sz == 0)
					continue;

				typeOfElem = GetTypeOfElement(nElem);

				if (typeOfElem <= 30) // parallelepiped
				{
				//initialization of an object of class T_Brick
				for (i=0; i<8; i++)
				{
					t = GetGlobalVertex(nElem, i);
					x[i] = xyz[t][0];
					y[i] = xyz[t][1];
					z[i] = xyz[t][2];
				}
				T_Brick L(x, y, z, typeOfElem);

				// determine the weights for the basis functions
				for (i=0; i<12; i++)
					ves[i] = v3[ed[nElem][i]];

				// for each point that fell into the element, we give the required value from the parallelepiped
				for (i=0; i<sz; i++)
				{
					t = PointresForElem[nElem][i];
					tx = pointres[t][0];
					ty = pointres[t][1];
					tz = pointres[t][2];

						switch(derive)
						{
						case 0: // Ax
							L.VectorFieldOnPar(tx,ty,tz, ves, &result[t], &tmp1, &tmp2);
							break;
						case 1: // Ay
							L.VectorFieldOnPar(tx,ty,tz, ves, &tmp1, &result[t], &tmp2);
							break;
						case 2: // Az
							L.VectorFieldOnPar(tx,ty,tz, ves, &tmp1, &tmp2, &result[t]);
							break;
						case 3: // dBz
							L.RotZOnPar(tz, ves, &result[t]);
							break;
						case 4: // dBx
							L.RotXOnPar(tx, ves, &result[t]);
							break;
						case 5: // dBy
							L.RotYOnPar(ty, ves, &result[t]);
							break;
						}
					}
				} 
					else // hexahedron
					{
						}
						}
					}
						}

	return 0;
	}
//------------------------------------------------------------------------
long Output3dArbitrary::GetGlobalVertex(long nElem, long nLocVertex)
{
	if (nvtr!=NULL)
	{
		return nvtr[nElem][nLocVertex];
	}
	else
	{
		return nver[nElem][nLocVertex];
	}
}
//------------------------------------------------------------------------
long Output3dArbitrary::GetTypeOfElement(long nElem)
{
	if (type_of_hex!=NULL)
	{
		return type_of_hex[nElem];
	}
	else
	{
		return nver[nElem][13];
	}
}
//------------------------------------------------------------------------
// 
//------------------------------------------------------------------------
void Output3dArbitrary::PointresSort()
{
	vector<PointRes2> p;
	long i;

	p.resize(n_pointres);

	// x
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
	
	// y
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

	// z
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
}
//------------------------------------------------------------------------
int Output3dArbitrary::IsPointInsideBrick(double xm, double ym, double zm, double *x0, double *x1)
{
	if (xm >= x0[0] && ym >= x0[1] && zm >= x0[2] &&
		xm <= x1[0] && ym <= x1[1] && zm <= x1[2])
	{
		return 1;
	}
	return 0;
}
//------------------------------------------------------------------------
int Output3dArbitrary::FindElemForReceivers()
{
	long i, j, sz, t;
	long typeOfHex;

	// boundaries of arrays with receivers, which fall (by x, y, z, respectively) into a parallelepiped outlined around the element
	vector<long_double>::iterator beg_x, beg_y, beg_z, end_x, end_y, end_z; 

	// max / min coordinates of the vertex of the element by x,y,z
	long_double coordMin[3], coordMax[3]; 

	vector<long_double> pntInBrick, pntInBrickTmp;
	vector<long_double> tmp_x, tmp_y, tmp_z;


	PointresForElem.resize(kpar);

	elemForPoint.resize(n_pointres);
	for(i=0; i<n_pointres; i++)
		elemForPoint[i] = -1;

	// for each element we find all the receivers that fall into it
	for(i=0; i<kpar; i++)
	{
		typeOfHex = GetTypeOfElement(i);

		if(typeOfHex <= 30) // parallelepiped
		{
			for (j=0; j<3; j++)
			{
				coordMin[j].d = xyz[GetGlobalVertex(i,0)][j];
				coordMax[j].d = xyz[GetGlobalVertex(i,7)][j];
			}

			// beg
			beg_x = lower_bound(PointresXsorted.begin(), PointresXsorted.end(), coordMin[0], Long_double_less());
			if(beg_x == PointresXsorted.end())
				continue;

			beg_y = lower_bound(PointresYsorted.begin(), PointresYsorted.end(), coordMin[1], Long_double_less());
			if(beg_y == PointresYsorted.end())
				continue;

			beg_z = lower_bound(PointresZsorted.begin(), PointresZsorted.end(), coordMin[2], Long_double_less());
			if(beg_z == PointresZsorted.end())
				continue;

			// end
			end_x = upper_bound(PointresXsorted.begin(), PointresXsorted.end(), coordMax[0], Long_double_less());
			if(end_x == PointresXsorted.begin())
				continue;

			end_y = upper_bound(PointresYsorted.begin(), PointresYsorted.end(), coordMax[1], Long_double_less());
			if(end_y == PointresYsorted.begin())
				continue;

			end_z = upper_bound(PointresZsorted.begin(), PointresZsorted.end(), coordMax[2], Long_double_less());
			if(end_z == PointresZsorted.begin())
				continue;

			// did not overlap beg and end
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

				if (elemForPoint[t] == -1)
				{
					elemForPoint[t] = i;
					PointresForElem[i].push_back(t);
				}
			}

			pntInBrick.clear();
			pntInBrickTmp.clear();
		}
	}

	// checking if there are receivers that did not hit any element
	for (i=0; i<n_pointres; i++)
		{
		if (elemForPoint[i]==-1)
			{
				char s[256];
			sprintf(s,"There is no element for receiver N %d",i);
				logfile<<s<<endl;
		}
	}
	
	return 0;
}
//------------------------------------------------------------------------