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
 *  This file contains a structures of finite element mesh storing                                   
 *                                                                                                   
 *  Written by Prof. Marina G. Persova, Ph.D. Petr A. Domnikov and Ph.D. Dmitry S. Kiselev                                       
 *  Novosibirsk State Technical University,                                                          
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                
 *  mpersova@mail.ru, p_domnikov@mail.ru, harlequin_00@mail.ru                                                         
 *  Version 2.1 December 17, 2020                                                                    
*/                                                                                                   


#pragma once
#include "T_Mapping.h"

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

	int Read_prep_data(); // reading the grid for VFEM MT
	int ReadPrepDataHarmLoop(char *pointres_fname); // read grid with harmonic source
	int Read_mtz_1d(); //reads usin.dat, ucos.dat, alfa, nu (this is for MT)      

	int Read_mesh_for_nonstat_problem(char *pointres_fname);
	int Read_infite0(); 

	
	int Read_3dmeshregular(long interval);  // reads 3dmeshregular skipping grid coordinates at "interval" (if we want to build a spline with a step larger than in 3dmeshregular)

	int LoadReceivers(char *fname, int& _n_pointres, double (*(&_pointres))[3]);

	long n_pointres;
	double (*pointres)[3]; // receiver coordinates    

	long maxiter; // maximum number of iterations for the solver (from config)
	double eps;   // epsilon for solver (from config)
	long n_materials; 	// number of different materials                                  
	double *mu3d;     	//(mu for 3D problem with objects (first column from mu3d file))  
	double *mu0;      	//(mu for one-dimensional (second column from mu3d file))         
	int n_pointresB;  	// number of receivers for B                                      
	int n_pointresE;  	// number of receivers for E                                      
	double (*pointresB)[3]; // receiver coordinates                                           
	double (*pointresE)[3]; // receiver coordinates                                           
	double *sigma3d;       	 //(sigma for a 3D problem with objects (first column from sig3d))
	double *sigma0;        	 //(sigma for one-dimensional (second column from sig3d))         
	double *dpr3d;      	//(dpr for a 3D problem with objects (first column from dpr3d))   
	double *dpr0;        	 //(dpr for one-dimensional (second column from dpr3d))           
	long kuzlov;     	 // number of nodes (all, including terminal ones)                
	long kpar;        	 // number of elements in the grid                                
	long kt1;         	 // number of nodes with first edge                               
	long *l13d;       	 // numbers of nodes with the first edge                          
	long (*nver)[14]; 	 // finite element node numbers (13-node hexagons)                
	long *nvkat;       	 // numbers of materials of final elements                        
	double (*xyz)[3];  	 // node coordinates
                                              
	double (*xyzt)[3];  	// node coordinates   
	long n_layers_1d;  	// number of different layers (materials) for a one-dimensional problem (from sreda1d.ay)    
	double *layers_1d; 	// layers for 1D problem (from sreda1d.ay)                                                   
	double *sigma_1d;  	// sigma0 for 1D problem (from sreda1d.ay)                                                   

	int LoadVectorE0ForLine(int n);
	vector< vector<EnLine> > EnForLine;

	int N_X,N_Y,N_Z;
	vector<int> regular;
	vector<double> Xcrd,Ycrd,Zcrd;

	// 3dmeshregular
	long n_mesh_regular_x;  // number of steps by x in 3dmeshregular             
	long n_mesh_regular_y;  // number of steps in y in 3dmeshregular             
	double *mesh_regular_x; // x-coordinates from 3dmeshregular                  
	double *mesh_regular_y; // y-coordinates from 3dmeshregular                  

	
	long ntime;     
	double *time;   

	
 	
 	double nu;    // frequency                                                                      
 	long alfa;    // current direction: by x - alpha=1; by y - alpha=0; J=(alpha*Jx,(1-alpha)*Jy, 0)
 	double *usin; // solution of a one-dimensional problem (sin-component)                          
 	double *ucos; // solution of a one-dimensional problem (sin-component)                          
 	double *z_1d; // one-dimensional grid                                                           
 	long n_1d;    // number of nodes in a one-dimensional grid                                      

	T_Mapping_Vec *tmap;
	long tasktype;

	int npr, nfreq;

	int fdirect;
};
//-----------------------------------------------------------