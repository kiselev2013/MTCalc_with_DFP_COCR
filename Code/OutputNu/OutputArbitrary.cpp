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
// дифференцирование решения по времени (по 3-слойной схеме) в точке t из [t_j2, t_j]
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
// конструктор для векторной задачи
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
// конструктор для узловой задачи
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
// выдать значения из сетки в точках и продифференцировать по 3-слойной схеме
// Для выдачи E=-dA/at и Eds=-dBz/dt
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
	double x[8], y[8], z[8]; // координаты текущего шестигранника, с к-рого выдаём рез-тат
	double tx, ty, tz; // координаты текущего приемника
	double ves_j[12], ves_j1[12], ves_j2[12]; // веса для ребер текущго шестигранника на трёх временных слоях
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
				sz = (long)PointresForElem[nElem].size(); // число приёмников, попавших в эл-т
				if(sz == 0)
					continue;

				typeOfElem = GetTypeOfElement(nElem);

				if (typeOfElem <= 30) // параллелепипед
				{
					// инициализация объекта класса T_Brick
					for (i=0; i<8; i++)
					{
						j = GetGlobalVertex(nElem, i);
						x[i] = xyz[j][0];
						y[i] = xyz[j][1];
						z[i] = xyz[j][2];
					}
					T_Brick L(x, y, z, typeOfElem);

					// определяем веса для баз. функций
					for (i=0; i<12; i++)
					{
						ves_j2[i] = v3_j2[ed[nElem][i]];
						ves_j1[i] = v3_j1[ed[nElem][i]];
						ves_j[i]  = v3_j[ed[nElem][i]];
					}

					// для каждой точки, попавшей в эл-т, выдаём с параллелепипеда требуемую величину
					for (i=0; i<sz; i++)
					{
						j = PointresForElem[nElem][i];
						tx = pointres[j][0];
						ty = pointres[j][1];
						tz = pointres[j][2];

						// выдаём значения на трёх временных слоях
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
						
						// собственно, дифференцируем по времени
						result[j] = -dA_dt(t, val_j, val_j1, val_j2, dt, dt0, dt1, t_j, t_j1, t_j2);
					}
				} 
				else // шестигранник
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
				sz = (long)PointresForElem[nElem].size(); // число приёмников, попавших в эл-т
				if(sz == 0)
					continue;

				typeOfElem = GetTypeOfElement(nElem);

				if (typeOfElem <= 30) // параллелепипед
				{
				// инициализация объекта класса T_Brick
				for (i=0; i<8; i++)
				{
					t = GetGlobalVertex(nElem, i);
					x[i] = xyz[t][0];
					y[i] = xyz[t][1];
					z[i] = xyz[t][2];
				}
				T_Brick L(x, y, z, typeOfElem);

				// определяем веса для баз. функций
				for (i=0; i<12; i++)
					ves[i] = v3[ed[nElem][i]];

				// для каждой точки, попавшей в эл-т, выдаём с параллелепипеда требуемую величину
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
					else // шестигранник
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

	// границы массивов с приёмниками, к-рые попадают (по x,y,z соответственно)
	// в параллелепипед, очерченный вокруг эл-та
	vector<long_double>::iterator beg_x, beg_y, beg_z, end_x, end_y, end_z; 

	// координаты макс/мин вершины эл-та по x,y,z
	long_double coordMin[3], coordMax[3]; 

	vector<long_double> pntInBrick, pntInBrickTmp;
	vector<long_double> tmp_x, tmp_y, tmp_z;


	PointresForElem.resize(kpar);

	elemForPoint.resize(n_pointres);
	for(i=0; i<n_pointres; i++)
		elemForPoint[i] = -1;

	// для каждого эл-та находим все приёмники, к-рые в него попадают
	for(i=0; i<kpar; i++)
	{
		typeOfHex = GetTypeOfElement(i);

		if(typeOfHex <= 30) // параллелепипед
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

			// не перехлестнулись ли beg и end
			if(distance(beg_x, end_x) <= 0)
				continue;

			if(distance(beg_y, end_y) <= 0)
				continue;

			if(distance(beg_z, end_z) <= 0)
				continue;

			// сейчас точки отсортированы по координатам
			// то что находится между beg_? и end_? необходимо отсортировать по номерам приёмников 
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

	// проверка есть ли приёмники, к-рые не попали ни в один эл-т
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