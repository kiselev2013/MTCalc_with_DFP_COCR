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
 *  Written by Ph.D. Petr A. Domnikov                                                     
 *  Novosibirsk State Technical University,                                               
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                     
 *  p_domnikov@mail.ru                                                                    
 *  Version 2.1 December 17, 2020                                                             
*/                                                                                        

#pragma once
class Vec_Prep_Data
{
public:

	int n_materials;  // number of different materials
	double* mu3d; 	 //(mu for 3D problem with objects (first column from mu3d file))
	double* mu0;	 //(mu for one-dimensional (second column from mu3d file))
	double* sigma3d; //(sigma for a 3D problem with objects (first column from sig3d))
	double* sigma0;	 //(sigma for one-dimensional (second column from sig3d))
	double* dpr3d;	//(dpr for a 3D problem with objects (first column from dpr3d))
	double* dpr0;	 //(dpr for one-dimensional (second column from dpr3d))
	int kuzlov;	 // number of nodes (all, including terminal ones)
	int kpar;	 // number of elements in the grid
	int kt1;	 // number of nodes with first edge
	int* l13d;	 // numbers of nodes with the first edge
	int(*nver)[14];	 // finite element node numbers (13-node hexagons)
	int* nvkat;	 // numbers of materials of final elements
	double(*xyz)[3]; // node coordinates

	Vec_Prep_Data();
	~Vec_Prep_Data();

	int Read_prep_data(); // reading the grid for VFEM MT
};
//-----------------------------------------------------------