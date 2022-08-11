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


#pragma once
#include "T_Brick.h"

class VecBlockSLAE
{
public:
	long *ig, *jg;

	double *di_block; // diagonal blocks
	double *gg_block; // off-diagonal blocks
	long *idi; // storage addresses of the beginnings of diagonal blocks
	long *ijg; // storage addresses of beginnings of off-diagonal blocks

	double *pr; // global right-hand side vector

	long n_elem;    // number of elements in the mesh
	long n_edges;   // number of edges (total)
	long n;         // SLAE dimension
	long nb;        // number of blocks (block dimension of SLAE)
	long ig_n_1;    // ig(n+1)-1    

	long (*nver)[14];
	long (*ed)[25];
	double (*xyz)[3];
	long (*edges)[2];

	long n_of_materials;
	long *nvkat;
	double *dpr3d;
	double *dpr0;
	double *sigma3d;
	double *sigma0;
	double *mu3d;
	double *mu0;
	double omega;

	long alpha;
	long n_1d;
	double *z_1d;
	double *sin_1d;
	double *cos_1d;

	VecBlockSLAE(long *ig, long *jg, long *idi, long *ijg,long n_elem, long n_edges,	
		double (*xyz)[3], long (*nver)[14], long (*ed)[25], long (*edges)[2], 
		long *nvkat, double nu, double *mu3d, double *mu0, double *sigma3d, double *sigma0, double *dpr3d, double *dpr0,
		long alpha, long n_1d, double *z_1d, double *sin_1d, double *cos_1d);

	~VecBlockSLAE();

	void AsmBlockSLAE(Vec_Prep_Data *d);

	int WriteBlockSLAEtoFiles();

	void Add_to_di_block(T_Brick *LocElem, int i, int j, long nBlock, double mult);
	void Add_to_gg_block(T_Brick *LocElem, int i, int j, long nBlock, double mult);
	void Add_to_pr_block(T_Brick *LocElem, int i, long nBlock, double mult);
};


