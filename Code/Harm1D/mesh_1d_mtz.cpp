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
 *  1D mesh generator     
 *                                                                                         
 *  Written by Ph.D. Petr A. Domnikov                                                      
 *  Novosibirsk State Technical University,                                                
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                      
 *  p_domnikov@mail.ru                                                                     
 *  Version 1.2 October 23, 2020                                                          
*/                                                                                         


#include "stdafx.h"
#include "mesh_1d_mtz.h"
#include "in_out.h"
//-----------------------------------------------------------
Mesh_1d_mtz::Mesh_1d_mtz()
{
	coords = NULL;
	nvkat = NULL;
	layers = NULL;
	mu = NULL;
	sigma = NULL;
}
//-----------------------------------------------------------
Mesh_1d_mtz::~Mesh_1d_mtz()
{
	if(coords) {delete [] coords; coords=NULL;}
	if(nvkat) {delete [] nvkat; nvkat=NULL;}
	if(layers) {delete [] layers; layers=NULL;}
	if(mu) {delete [] mu; mu=NULL;}
	if(sigma) {delete [] sigma; sigma=NULL;};
}
//-----------------------------------------------------------
int Mesh_1d_mtz::Read_1d_data_for_3d_problem()
{
	FILE *fp;
	long i;
	double tmp1, tmp2;
	In_Out r;

	//sreda1d.ay
	if((fp=fopen("sreda1d.ay","r"))==0)
	{
		printf("Cannot open sreda1d.ay\n");
		return 1;
	}
	fscanf(fp, "%ld", &n_layers);

	n_materials = n_layers + 1 ; // air

	if((layers = new double[n_layers+1])==0) Memory_allocation_error("layers", "Mesh_1d_mtz::Read_1d_data_for_3d_problem");

	for(i=0; i<n_layers; i++)
		fscanf(fp, "%lf %lf %lf", &layers[i+1], &tmp1, &tmp2);

	fclose(fp);
	layers[0] = X0_1D;

	//mu3d
	if((mu = new double[n_materials])==0) Memory_allocation_error("mu", "Mesh_1d_mtz::Read_1d_data_for_3d_problem");

	for(i=0; i<n_materials; i++)
		mu[i] = 4.0*PI*1E-7;

	//sig3d
	if((sigma = new double[n_materials])==0) Memory_allocation_error("sigma", "Mesh_1d_mtz::Read_1d_data_for_3d_problem");;

	if((fp=fopen("sig3d", "r"))==0)
	{
		printf("Error: Cannot open file sig3d.\n");
		return 1;
	}

	for(i=0; i<n_materials; i++)
		fscanf(fp, "%lf %lf %lf", &tmp1, &tmp2, &sigma[i]);
	
	fclose(fp);

	//nu
	r.Read_Double_From_Txt_File("nu", &omega);
	omega *= 2.0*PI;

	return 0;  
}
//-----------------------------------------------------------
void Mesh_1d_mtz::Gen_1d_mesh()
{
	double x0;
	double x;
	double step;
	long i, j, k;
	double *temp;
	double center;

	x = x0 = X1_1D;
	step = STEP0;
	n_points = 0; 

	if((nodes_in_layer = new long[n_layers])==0) Memory_allocation_error("nodes_in_layer", "Mesh_1d_mtz::Gen_1d_mesh");

	for(i=0; i<n_layers; i++)
		nodes_in_layer[i] = 0;

	// first determine how many nodes will be in each layer
	for(i=0; i<n_layers; i++)
	{
		while(true)
		{
			x -= step;
			if(x <= layers[n_layers-i-1])
				break;
			step *=COEF_RAZR;
			nodes_in_layer[i]++;	
		}
	}

	// total number of nodes
	for(i=0; i<n_layers; i++)
		n_points += nodes_in_layer[i];
	n_points += n_layers + 1; // add nodes at layer interfaces

	n_elem = n_points - 1;
	if((coords = new double[n_points])==0) Memory_allocation_error("coords", "Mesh_1d_mtz::Gen_1d_mesh");

	if((nvkat = new long[n_elem])==0) Memory_allocation_error("nvkat", "Mesh_1d_mtz::Gen_1d_mesh");

	// determine the coordinates of the finite element nodes
	x = x0 = X1_1D;
	step = STEP0;
	k = 0;
	for(i=0; i<n_layers; i++)
	{
		coords[k] = layers[n_layers - i];
		k++;

		while(true)
		{
			x -= step;
			if(x <= layers[n_layers-i-1])
				break;
			step *=COEF_RAZR;
			coords[k] = x;	
			k++;
		}
	}
	coords[k] = X0_1D;

	if((temp = new double[n_points])==0) Memory_allocation_error("temp", "Mesh_1d_mtz::Gen_1d_mesh");

	for(i=0; i<n_points; i++)
		temp[i] = coords[n_points-i-1];

	for(i=0; i<n_points; i++)
		coords[i] = temp[i];

	delete [] temp;

	// nvkat
	for(i=0; i<n_elem; i++)
	{
		center = 0.5*(coords[i] + coords[i+1]);

		for(j=0; j<n_layers; j++)
		{
			if(center>layers[j] && center<layers[j+1])
			{
				nvkat[i] = j+1;
				break;
			}
		}
	}

	delete [] nodes_in_layer;

	In_Out r;
	r.Write_Txt_File_Of_Double("mesh1d", coords, n_points, 1);
	r.Write_Txt_File_Of_Long("nvkat1d", nvkat, n_elem, 1);
}
//-----------------------------------------------------------