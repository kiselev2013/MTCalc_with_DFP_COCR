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
 *  This file contains the routine to solve the 1D MT problem in the layered medium by FEM          
 *                                                                                                     
 *  Written by Ph.D. Petr A. Domnikov                                                                  
 *  Novosibirsk State Technical University,                                                            
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                  
 *  p_domnikov@mail.ru                                                                                 
 *  Version 1.3 November 23, 2020                                                                       
*/                                                                                                     

#include "stdafx.h"
#include "divgrad_1d.h"
#include "mesh_1d_mtz.h"
#include "portrait_harm_prof.h"
#include "solver_profil.h"
#include "global_slae_1d_harm_prof.h"
#include "in_out.h"
extern ofstream logfile;
//-----------------------------------------------------------
Divgrad_1d::Divgrad_1d()
{
	this->n_1d = 0;
	this->alpha = alpha;

	this->coords_1d = NULL;
	this->cos_1d = NULL;
	this->sin_1d = NULL;
}
//--------------------------------------------------------------
Divgrad_1d::~Divgrad_1d()
{
	if(cos_1d) {delete [] cos_1d; cos_1d=NULL;}
	if(sin_1d) {delete [] sin_1d; sin_1d=NULL;}
	if(coords_1d) {delete [] coords_1d; coords_1d=NULL;}
}
//-----------------------------------------------------------
struct Point3D
{
	double x,y,z;
};
double Spline(double x, long n, double *xyz, double *values)
{
	double s, xi;
	long i, t, flag;

	flag = 0;

	if (x < xyz[0])
	{
		return values[0];
	}

	if (x > xyz[n-1])
	{
		return values[n-1];
	}

	for(i=0; i<n-1; i++)
	{
		if(x >= xyz[i]  &&  x <= xyz[i+1])
		{
			t = i;
			flag = 1;
			break;
		}
	}

	if(flag == 1)
	{
		xi = (x - xyz[t])/(xyz[t+1] - xyz[t]);
		s = (1.0 - xi)*values[t] + xi*values[t+1];
	}
	else
	{
		s = 0.0;
	}

	return s;
}
//-----------------------------------------------------------
double dSpline(double x, long n, double *xyz, double *values)
{
	int i,fl[2];
	double coeff,df[2],ddf[2];
	if (x < xyz[1])
	{
		return (values[1]-values[0])/(xyz[1] - xyz[0]);
	}
	else if (x > xyz[n-2])
	{
		return (values[n-1]-values[n-2])/(xyz[n-1] - xyz[n-2]);
	}
	else
	{
		fl[0]=fl[1]=0.0;
		df[0]=df[1]=0.0;
		for(i=1;i<n; i++)
		{
			if(x <= xyz[i])
			{
				break;
			}
		}
		if(i-1>0)
		{
			fl[0]=1;
			coeff=(xyz[i-1]-xyz[i-2])/(xyz[i] - xyz[i-2]);
			ddf[0]=(values[i-1]-values[i-2])/(xyz[i-1] - xyz[i-2]);
			ddf[1]=(values[i]-values[i-1])/(xyz[i] - xyz[i-1]);
			df[0]=ddf[0]*(1.0-coeff)+ddf[1]*coeff;
		}
		if(i+1<n)
		{
			fl[1]=1;
			coeff=(xyz[i]-xyz[i-1])/(xyz[i+1] - xyz[i-1]);
			ddf[0]=(values[i]-values[i-1])/(xyz[i] - xyz[i-1]);
			ddf[1]=(values[i+1]-values[i])/(xyz[i+1] - xyz[i]);
			df[1]=ddf[0]*(1.0-coeff)+ddf[1]*coeff;
		}
		if(fl[0] && fl[1])
		{
			coeff=(x-xyz[i-1])/(xyz[i] - xyz[i-1]);
			return (1.0-coeff)*df[0]+coeff*df[1];
		}
		else
		{
			return (values[i]-values[i-1])/(xyz[i] - xyz[i-1]);
		}
	}
}
//-----------------------------------------------------------
ofstream& operator<(ofstream& file, const double &id)
{
	file.write((char*)&id, sizeof(double));
	return file;
}
//-----------------------------------------------------------
int Divgrad_1d::Solve_1d_Problem_for_3d_task()
{
	Mesh_1d_mtz *m=NULL;
	Portrait_profil *p=NULL;
	Global_slae_1d_harm_prof *a=NULL;
	Solver_profil *s=NULL;
	long n;
	long i;
	double *solution=NULL;
	FILE *fp=NULL;
	In_Out r;

	logfile<<1<<endl;

	// generate a one-dimensional grid
	m = new Mesh_1d_mtz();
	if(m == 0)
		Memory_allocation_error("m", "Divgrad_1d::Solve_1d_Problem_for_3d_task");

	m->Read_1d_data_for_3d_problem();
	m->Gen_1d_mesh();

	// matrix portrait
	p = new Portrait_profil(m->n_elem);
	if(p == 0)
		Memory_allocation_error("p", "Divgrad_1d::Solve_1d_Problem_for_3d_task");

	p->Gen_ig_for_1d_harm_problem();

	n = m->n_points*2;

	// SLAE assembly
	a = new Global_slae_1d_harm_prof(n, m->n_elem, p->ig, m->nvkat, m->mu,
		m->sigma, m->omega, m->coords, 0);
	if(a == 0)
		Memory_allocation_error("a", "Divgrad_1d::Solve_1d_Problem_for_3d_task");

	a->Assembling_for_1d_harm_problem();

	// solve SLAE
	solution = new double[n];
	s = new Solver_profil(n, p->ig, a->di, a->ggl, a->ggu, a->pr, solution); 
	if(s == 0)
		Memory_allocation_error("s", "Divgrad_1d::Solve_1d_Problem_for_3d_task");

	s->Solve_SLAE_using_LU();

	// writing the result
	// usin.dat
	if((fp=fopen("usin.dat", "w"))==0)
		Cannot_open_file("usin.dat", "Divgrad_1d::Solve_1d_Problem_for_3d_task");

	for(i=0; i<m->n_points; i++)
		fprintf(fp, "%25.13e %25.13e\n", m->coords[i], solution[i*2]);

	fclose(fp);

	// ucos.dat
	if((fp=fopen("ucos.dat", "w"))==0)
		Cannot_open_file("ucos.dat", "Divgrad_1d::Solve_1d_Problem_for_3d_task");

	for(i=0; i<m->n_points; i++)
		fprintf(fp, "%25.13e %25.13e\n", m->coords[i], solution[i*2+1]);

	fclose(fp);

	sin_1d=new double[m->n_points];
	cos_1d=new double[m->n_points];

	// setka1DEy
	r.Write_Long_To_Txt_File("setka1DEy", m->n_points);
	
	r.Read_1d_data(m->n_points,m->coords,sin_1d,cos_1d);


	ifstream inf;
	ofstream ofp;
	int npntE,npntE0,npntB;
	vector<Point3D> pntE,pntB;

	npntE=npntE0=npntB=0;
	pntE.clear();
	pntB.clear();

	inf.open("xyzVectorB");
	if(inf)
	{
		inf>>npntB;
		pntB.resize(npntB);
		for(i=0;i<npntB;i++)
		{
			inf>>pntB[i].x>>pntB[i].y>>pntB[i].z;
		}
		inf.close();
	}
	inf.clear();

	inf.open("xyzVectorE");
	if(inf)
	{
		inf>>npntE;
		pntE.resize(npntE);
		for(i=0;i<npntE;i++)
		{
			inf>>pntE[i].x>>pntE[i].y>>pntE[i].z;
		}
		inf.close();
	}
	inf.clear();

	inf.open("xyzVectorE0");
	if(inf)
	{
		inf>>npntE0;
		pntE.resize(npntE+npntE0);
		for(i=0;i<npntE0;i++)
		{
			inf>>pntE[npntE+i].x>>pntE[npntE+i].y>>pntE[npntE+i].z;
		}
		inf.close();
	}
	inf.clear();


	ofstream outfs,outfc;

	int alfa;
	double alpha_double;
	double usin, ucos, ducos, dusin;
	double w, mu;

	r.Read_Double_From_Txt_File("alfa", &alpha_double);
	if (fabs(1.0 - alpha_double)<0.01)
	{
		alfa = 1; 
	} 
	else
	{
		alfa = 0;
	}

	w = m->omega;
	mu = 4.0*PI*1E-7;


	ofp.open("b2d");
	ofp<<scientific<<setprecision(14);
	for(i=0;i<npntB;i++)
	{
		if(pntB[i].z<0.0)
		{
			dusin = dSpline(pntB[i].z, m->n_points, m->coords, sin_1d);
			ducos = dSpline(pntB[i].z, m->n_points, m->coords, cos_1d);
		}
		else
		{
			dusin = mu;
			ducos = 0.0;
		}
		ofp<<(alfa - 1.0)*dusin<<" "<<alfa*dusin<<" "<<0.0<<" "<<
			(alfa - 1.0)*ducos<<" "<<alfa*ducos<<" "<<0.0<<'\n';
	}
	ofp.close();
	ofp.clear();


	ofp.open("e2d");
	ofp<<scientific<<setprecision(14);
	for(i=0;i<npntE;i++)
	{
		if(!(pntE[i].z>0.0))
		{
			usin = Spline(pntE[i].z, m->n_points, m->coords, sin_1d);
			ucos = Spline(pntE[i].z, m->n_points, m->coords, cos_1d);
			ofp<<w*alfa*ucos<<" "<<w*(1.0 - alfa)*ucos<<" "<<0.0<<" "<<
				-w*alfa*usin<<" "<<w*(alfa - 1.0)*usin<<" "<<0.0<<'\n';
		}
		else
		{
			//ofp<<0.0<<" "<<0.0<<" "<<0.0<<" "<<
			//	0.0<<" "<<0.0<<" "<<0.0<<'\n';
			usin = Spline(0.0, m->n_points, m->coords, sin_1d);
			ucos = Spline(0.0, m->n_points, m->coords, cos_1d);
			ofp<<w*alfa*ucos<<" "<<w*(1.0 - alfa)*ucos<<" "<<0.0<<" "<<
				-w*alfa*usin<<" "<<w*(alfa - 1.0)*usin<<" "<<0.0<<'\n';
		}
	}
	ofp.close();
	ofp.clear();


	if(npntE0)
	{
		outfs.open("e_s.dat", ios::binary);
		outfc.open("e_c.dat", ios::binary);
		for(i=0;i<npntE0;i++)
		{
			if(!(pntE[i].z>0.0))
			{
				usin = Spline(pntE[npntE+i].z, m->n_points, m->coords, sin_1d);
				ucos = Spline(pntE[npntE+i].z, m->n_points, m->coords, cos_1d);
				outfs<w*alfa*ucos<w*(1.0 - alfa)*ucos<0.0;
				outfc<-w*alfa*usin<w*(alfa - 1.0)*usin<0.0;
			}
			else
			{
				//outfs<0.0<0.0<0.0;
				//outfc<0.0<0.0<0.0;
				usin = Spline(0.0, m->n_points, m->coords, sin_1d);
				ucos = Spline(0.0, m->n_points, m->coords, cos_1d);
				outfs<w*alfa*ucos<w*(1.0 - alfa)*ucos<0.0;
				outfc<-w*alfa*usin<w*(alfa - 1.0)*usin<0.0;
			}
		}
		outfs.close();
		outfs.clear();
		outfc.close();
		outfc.clear();
	}


	if(m) {delete m; m=NULL;}
	if(p) {delete p; p=NULL;}
	if(a) {delete a; a=NULL;}
	if(s) {delete s; s=NULL;}
	if(solution) {delete [] solution; solution=NULL;}
	if(sin_1d) {delete [] sin_1d; sin_1d=NULL;}
	if(cos_1d) {delete [] cos_1d; cos_1d=NULL;}

	return 0;
}
//-----------------------------------------------------------
