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
 *  This file contains some basic routines for assembling the finite element matrix and vector of the right hand side
 *                                                                                                                   
 *  Written by Ph.D. Petr A. Domnikov                                                                                
 *  Novosibirsk State Technical University,                                                                          
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                                
 *  p_domnikov@mail.ru                                                                                               
 *  Version 1.3 April 5, 2021                                                                                        
*/                                                                                                                   

#include "stdafx.h"
#include "t_global_slae.h"

extern ofstream logfile;

VecBlockSLAE::VecBlockSLAE(long *ig, long *jg, long *idi, long *ijg, 
			  long n_elem, long n_edges,	
			  double (*xyz)[3], long (*nver)[14], long (*ed)[25], long (*edges)[2], 
			  long *nvkat, double nu, double *mu3d, double *mu0, double *sigma3d, double *sigma0, double *dpr3d, double *dpr0,
			  long alpha, long n_1d, double *z_1d, double *sin_1d, double *cos_1d)
{
	this->ig = ig;
	this->jg = jg;
	this->idi = idi;
	this->ijg = ijg;

	this->n_elem = n_elem;
	this->n_edges = n_edges; 

	this->xyz = xyz;
	this->nver = nver;
	this->ed = ed;
	this->edges = edges;

	this->nvkat = nvkat;
	this->omega = 2.0*PI*nu;
	this->mu3d = mu3d;
	this->mu0 = mu0;
	this->sigma3d = sigma3d;
	this->sigma0 = sigma0;
	this->dpr3d = dpr3d;
	this->dpr0 = dpr0;

	this->alpha = alpha;
	this->n_1d = n_1d;
	this->z_1d = z_1d;
	this->sin_1d = sin_1d;
	this->cos_1d = cos_1d;

	this->nb = n_edges;
	this->n = n_edges*2;

	this->ig_n_1 = this->ig[this->n_edges];

	di_block = NULL;
	gg_block = NULL;
	pr = NULL;

	if((pr = new double[this->n])==0) Memory_allocation_error("pr", "VecBlockSLAE::VecBlockSLAE");
	if((di_block = new double[idi[n_edges]])==0) Memory_allocation_error("di_block", "VecBlockSLAE::VecBlockSLAE");
	if((gg_block = new double[ijg[ig_n_1]])==0) Memory_allocation_error("gg_block", "VecBlockSLAE::VecBlockSLAE");
}

VecBlockSLAE::~VecBlockSLAE()
{
	if(pr)  { delete [] pr; pr=NULL; }
	if(gg_block) {delete [] gg_block; gg_block=NULL;}
	if(di_block) {delete [] di_block; gg_block=NULL;}
}

void VecBlockSLAE::AsmBlockSLAE(Vec_Prep_Data *d)
{
	long i, j, k, m, it, jt, i_mu, j_nu;
	long ii, jj; // global numbers (for T-transform matrix)

	//set 0
	for(i=0; i<=idi[nb]-1; i++)
		di_block[i] = 0;

	for(i=0; i<=ijg[ig_n_1]-1; i++)
		gg_block[i] = 0;

	for(i=0; i<n; i++)
		pr[i] = 0.0;

	for(i=0; i<n_elem; i++)
	{
		T_Brick L(i, nver, ed, edges, xyz, nvkat, sigma3d, sigma0, mu3d, omega,
			alpha, n_1d, z_1d, sin_1d, cos_1d);

		L.d=d;

		for(j=0;j<24;j++){L.g_harm[j]=0.0;}

		L.mu0=mu0[nvkat[i]];
		L.mu=mu3d[nvkat[i]];			
		L.sigma0=sigma0[nvkat[i]];
		L.sigma=sigma3d[nvkat[i]];
		L.dpr0=dpr0[nvkat[i]];
		L.dpr=dpr3d[nvkat[i]];

		L.Calc_block_local_matrix_and_vector();

		for(j=0; j<12; j++)
		{
			ii = ed[i][j];

			Add_to_di_block(&L, j,j, ii, 1.0); 
			Add_to_pr_block(&L, j, ii, 1.0);

			for(k=0; k<12; k++)
			{
				jj = ed[i][k];

				if(jj < ii) 
					for(m=ig[ii]; m<=ig[ii+1]-1; m++)
						if(jg[m]==jj)
						{
							Add_to_gg_block(&L, j,k, m, 1.0);
							break;
						}
			}//k
		}//j
	}//i
}

void VecBlockSLAE::Add_to_di_block(T_Brick *LocElem, int i, int j, long nBlock, double mult)
{
	long beg = idi[nBlock];
	long size = idi[nBlock+1] - beg;

	if(size==1)
	{
		di_block[beg] += LocElem->b[i][j]*mult;  	
	}
	else
	{
		di_block[beg] += LocElem->b[i][j]*mult;
		di_block[beg+1] += LocElem->c[i][j]*mult;  	
	}
}

void VecBlockSLAE::Add_to_gg_block(T_Brick *LocElem, int i, int j, long nBlock, double mult)
{
	long beg = ijg[nBlock];
	long size = ijg[nBlock+1] - beg;

	if(size==1)
	{
		gg_block[beg] += LocElem->b[i][j]*mult;  	
	}
	else
	{
		gg_block[beg] += LocElem->b[i][j]*mult;
		gg_block[beg+1] += LocElem->c[i][j]*mult;  	
	}
}

void VecBlockSLAE::Add_to_pr_block(T_Brick *LocElem, int i, long nBlock, double mult)
{
	pr[nBlock*2]   += LocElem->g_harm[i*2]*mult;
	pr[nBlock*2+1] += LocElem->g_harm[i*2+1]*mult;
}

int VecBlockSLAE::WriteBlockSLAEtoFiles()
{
	int maxiter=10000;
	double maxnev=1e-6;
	In_Out R;
	ifstream inf;

	inf.open("issc.txt");
	if(inf)
	{
		inf>>maxiter;
		inf>>maxnev;
		inf.close();
	}
	inf.clear();

	R.Write_kuslau("kuslau", n, maxnev, maxiter);

	if(pr) R.Write_Bin_File_Of_Double("pr", pr, n, 1);
	if(di_block) R.Write_Bin_File_Of_Double("di", di_block, idi[n_edges], 1);
	if(gg_block) R.Write_Bin_File_Of_Double("gg", gg_block, ijg[ig_n_1], 1);

	R.Write_Bin_File_Of_Long("ig", ig, n_edges+1, 1);
	R.Write_Bin_File_Of_Long("jg", jg, ig_n_1, 1);
	R.Write_Bin_File_Of_Long("idi", idi, n_edges+1, 1);
	R.Write_Bin_File_Of_Long("ijg", ijg, ig_n_1+1, 1);
	
	return 0;
}
