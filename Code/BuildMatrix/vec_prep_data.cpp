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
 *  Written by Ph.D. Petr A. Domnikov and Ph.D. Dmitry S. Kiselev                                                    
 *  Novosibirsk State Technical University,                                               
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                     
 *  Corresponding author:                                                                                            
 *  E-mail: p_domnikov@mail.ru (Petr A. Domnikov) 
 *  Version 2.1 December 17, 2020                                                         
*/                                                                                        
#include "stdafx.h" 
#include "vec_prep_data.h"

extern ofstream logfile;

Vec_Prep_Data::Vec_Prep_Data()
{
	nver = NULL;
	nvkat = NULL;
	xyz = NULL;
	mu3d = NULL;
	mu0 = NULL;
    sigma3d = NULL;
	sigma0 = NULL;
	n_pointresB=0;
	pointresB = NULL;
	n_pointresE=0;
	pointresE = NULL;
	mesh_regular_x = NULL;
	mesh_regular_y = NULL;
	l13d = NULL;
	usin = NULL;
	ucos = NULL;
	z_1d = NULL;
	sigma_1d = NULL;
	layers_1d = NULL;
	dpr3d = NULL;
	dpr0 = NULL;
}

Vec_Prep_Data::~Vec_Prep_Data()
{
	if(mu0) {delete [] mu0; mu0=NULL;}
	if(xyz) {delete [] xyz; xyz=NULL;}
	if(nver) {delete [] nver; nver=NULL;}
	if(mu3d) {delete [] mu3d; mu3d=NULL;}
	if(l13d) {delete [] l13d; l13d=NULL;}
	if(usin) {delete [] usin; usin=NULL;}
	if(ucos) {delete [] ucos; ucos=NULL;}
	if(z_1d) {delete [] z_1d; z_1d=NULL;}
	if(nvkat) {delete [] nvkat; nvkat=NULL;}
	if(sigma0) {delete [] sigma0; sigma0=NULL;}
	if(sigma3d) {delete [] sigma3d; sigma3d=NULL;}
	if(pointresB) {delete [] pointresB; pointresB=NULL;}
	if(pointresE) {delete [] pointresE; pointresE=NULL;}
	if(sigma_1d) {delete [] sigma_1d; sigma_1d=NULL;}
	if(layers_1d) {delete [] layers_1d; layers_1d=NULL;}
	if(mesh_regular_x) {delete [] mesh_regular_x; mesh_regular_x=NULL;}
	if(mesh_regular_y) {delete [] mesh_regular_y; mesh_regular_y=NULL;}
	if(dpr3d) {delete [] dpr3d; dpr3d=NULL;}
	if(dpr0) {delete [] dpr0; dpr0=NULL;}
}
//-----------------------------------------------------------------------------  
// reading the grid for VFEM MT                                                  
//-----------------------------------------------------------------------------  
int Vec_Prep_Data::Read_prep_data()
{
	In_Out R;
	FILE *fp=NULL;
	double temp, temp1, temp2;
	long i, j;
	long max_material;
	long tmp;
	double ang,sa,ca,ma[3][3],tm[3][3];

	int n_reg;
	ifstream inf;
	inf.open("regular",ios::binary);
	if (!inf)
		Cannot_open_file("regular", "Vec_Prep_Data::Read_prep_data");
	inf>n_reg;
	regular.resize(n_reg);
	for(i=0;i<n_reg;i++){
		inf>regular[i];
		regular[i]--;
	}
	inf.close();
	inf.clear();

	if((fp=fopen("sig3d", "r"))==0)
		Cannot_open_file("sig3d", "Vec_Prep_Data::Read_prep_data");

	max_material = 0;
	while(!feof(fp))
	{
		fscanf(fp, "%ld %lf %lf", &tmp, &temp1, &temp2);
		if(tmp > max_material)
			max_material = tmp;
	}
	n_materials = max_material;
	fclose(fp);

	mu3d = new double[n_materials];
	if(mu3d == 0)
		Memory_allocation_error("mu3d", "Vec_Prep_Data::Read_prep_data");

	mu0 = new double[n_materials];
	if(mu0 == 0)
		Memory_allocation_error("mu0", "Vec_Prep_Data::Read_prep_data");

	if((fp=fopen("mu3d", "r"))==0)
		Cannot_open_file("mu3d", "Vec_Prep_Data::Read_prep_data");

	for (i=0; i<n_materials; i++)
	{
		long n_of_current_material;
		fscanf(fp, "%ld", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf %lf",
			&mu3d[n_of_current_material], &mu0[n_of_current_material]);
		mu3d[n_of_current_material]*=MU_0;
		mu0[n_of_current_material]*=MU_0;
	}
	fclose(fp);

	sigma3d = new double[n_materials];
	if(sigma3d == 0)
		Memory_allocation_error("sigma3d", "Vec_Prep_Data::Read_prep_data");

	sigma0 = new double[n_materials];
	if(sigma0 == 0)
		Memory_allocation_error("sigma0", "Vec_Prep_Data::Read_prep_data");	

	if((fp=fopen("sig3d", "r"))==0)
		Cannot_open_file("sig3d", "Vec_Prep_Data::Read_prep_data");

	for (i=0; i<n_materials; i++)
	{
		long n_of_current_material;
		fscanf(fp, "%ld", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf %lf",
				&sigma3d[n_of_current_material], &sigma0[n_of_current_material]);
	}
	fclose(fp);

	dpr3d = new double[n_materials];
	if(dpr3d == 0)
		Memory_allocation_error("dpr3d", "Vec_Prep_Data::Read_prep_data");

	dpr0 = new double[n_materials];
	if(dpr0 == 0)
		Memory_allocation_error("dpr0", "Vec_Prep_Data::Read_prep_data");

	if((fp=fopen("dpr3d", "r"))==0)
		Cannot_open_file("dpr3d", "Vec_Prep_Data::Read_prep_data");

	for (i=0; i<n_materials; i++)
	{
		long n_of_current_material;
		fscanf(fp, "%ld", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf %lf",
			&dpr3d[n_of_current_material], &dpr0[n_of_current_material]);
		dpr3d[n_of_current_material]*=DPR_0;
		dpr0[n_of_current_material]*=DPR_0;
	}
	fclose(fp);

	LoadReceivers("pointres", n_pointresB, pointresB);
	LoadReceivers("pointres", n_pointresE, pointresE);

	R.Read_inftry("inftry.dat", &kuzlov, &kpar, &kt1);

	nver = new long[kpar][14];
	if(nver == 0)
		Memory_allocation_error("nver", "Vec_Prep_Data::Read_prep_data");

	R.Read_Bin_File_Of_Long("nver.dat", (long*)nver, kpar, 14);
	for(i=0; i<kpar; i++)
		for(j=0; j<14; j++)
			nver[i][j]--;

	nvkat = new long[kpar];
	if(nvkat == 0)
		Memory_allocation_error("nvkat", "Vec_Prep_Data::Read_prep_data");

	R.Read_Bin_File_Of_Long("nvkat.dat", nvkat, kpar, 1);
	for(i=0; i<kpar; i++)
		nvkat[i]--;

	xyz = new double[kuzlov][3];
	if(xyz == 0)
		Memory_allocation_error("xyz", "Vec_Prep_Data::Read_prep_data");

	R.Read_Bin_File_Of_Double("xyz.dat", (double*)xyz, kuzlov, 3);

	l13d = new long[kt1];
	if(l13d == 0)
		Memory_allocation_error("l13d", "Vec_Prep_Data::Read_prep_data");

	R.Read_Bin_File_Of_Long("l13d.dat", l13d, kt1, 1);
	for(i=0; i<kt1; i++)
		l13d[i]--;

	if(Read_3dmeshregular(0)!=0)
		throw logic_error("Error reading file 3dmeshregular.");

	if((fp=fopen("sreda1d.ay","r"))==0)
		Cannot_open_file("sreda1d.ay", "Vec_Prep_Data::Read_prep_data");

	fscanf(fp, "%ld", &n_layers_1d);

	layers_1d = new double[n_layers_1d];
	if(layers_1d == 0)
		Memory_allocation_error("layers_1d", "Vec_Prep_Data::Read_prep_data");

	sigma_1d = new double[n_layers_1d];
	if(sigma_1d == 0)
		Memory_allocation_error("sigma_1d", "Vec_Prep_Data::Read_prep_data");

	for(i=0; i<n_layers_1d; i++)
	{
		fscanf(fp, "%lf", &layers_1d[i]);
		fscanf(fp, "%lf", &sigma_1d[i]);
		fscanf(fp, "%lf", &temp);
	}

	fclose(fp);

	return 0;
}

int Vec_Prep_Data::Read_mtz_1d()
{
	In_Out r;
	double alpha_double;

	// alfa
	r.Read_Double_From_Txt_File("alfa", &alpha_double);
	if (fabs(1.0 - alpha_double)<0.01)
	{
		alfa = 1; 
	} 
	else
	{
		alfa = 0;
	}

	r.Read_Double_From_Txt_File("nu", &nu);
	r.Read_Long_From_Txt_File("setka1DEy", &n_1d);

	z_1d = new double[n_1d];
	if(z_1d == 0)
		Memory_allocation_error("z_1d", "Vec_Prep_Data::Read_mtz_1d");

	usin = new double[n_1d];
	if(usin == 0)
		Memory_allocation_error("usin", "Vec_Prep_Data::Read_mtz_1d");

	ucos = new double[n_1d];
	if(ucos == 0)
		Memory_allocation_error("ucos", "Vec_Prep_Data::Read_mtz_1d");

	return r.Read_1d_data(n_1d, z_1d, usin, ucos);
}

int Vec_Prep_Data::Read_3dmeshregular(long interval)
{
	FILE *fp=NULL;
	long i, j, t, r;
	double temp;
	long int_whole, int_ost;
	long nx, ny, nz;

	if((fp=fopen("3dmeshregular","r"))==0)
		Cannot_open_file("3dmeshregular", "Vec_Prep_Data::Read_3dmeshregular");

	fscanf(fp, "%ld", &nx);

	N_X=nx;
	Xcrd.resize(N_X);

	int_whole = nx/(1+interval);
	if(int_whole*(1+interval)-nx==0)
	{
		int_ost = 0;
	}
	else
	{	
		int_ost = 1;	
	}
	n_mesh_regular_x = int_whole + int_ost;

	mesh_regular_x = new double[n_mesh_regular_x];
	if(mesh_regular_x == 0)
		Memory_allocation_error("mesh_regular_x", "Vec_Prep_Data::Read_3dmeshregular");

	r=0;
	for(i=0; i<int_whole; i++)
	{
		t = i*(interval+1);
		fscanf(fp, "%lf %lf", &temp, &(Xcrd[r]));
		mesh_regular_x[i]=Xcrd[r];
		r++;
		for(j=0; j<interval; j++)
		{
			fscanf(fp, "%lf %lf", &temp, &(Xcrd[r]));
			r++;
			t++;
		}
	}

	if(int_ost!=0)
	{
		for(i=t+1; i<nx-1; i++){
			fscanf(fp, "%lf %lf", &temp, &(Xcrd[r]));
			r++;
		}
		fscanf(fp, "%lf %lf", &temp, &(Xcrd[r]));
		mesh_regular_x[n_mesh_regular_x-1]=Xcrd[r];
		r++;
	}

	fscanf(fp, "%ld", &ny);

	N_Y=ny;
	Ycrd.resize(N_Y);

	int_whole = ny/(1+interval);
	if(int_whole*(1+interval)-ny==0)
	{
		int_ost = 0;
	}
	else
	{	
		int_ost = 1;	
	}
	n_mesh_regular_y = int_whole + int_ost;

	mesh_regular_y = new double[n_mesh_regular_y];
	if(mesh_regular_y == 0)
		Memory_allocation_error("mesh_regular_y", "Vec_Prep_Data::Read_3dmeshregular");

	r=0;
	for(i=0; i<int_whole; i++)
	{
		t = i*(interval+1);
		fscanf(fp, "%lf %lf", &temp, &(Ycrd[r]));
		mesh_regular_y[i]=Ycrd[r];
		r++;
		for(j=0; j<interval; j++)
		{
			fscanf(fp, "%lf %lf", &temp, &(Ycrd[r]));
			r++;
			t++;
		}
	}

	if(int_ost!=0)
	{
		for(i=t+1; i<ny-1; i++){
			fscanf(fp, "%lf %lf", &temp, &(Ycrd[r]));
			r++;
		}
		fscanf(fp, "%lf %lf", &temp, &(Ycrd[r]));
		mesh_regular_y[n_mesh_regular_y-1]=Ycrd[r];
		r++;
	}

	fscanf(fp, "%ld", &nz);

	N_Z=nz;
	Zcrd.resize(N_Z);
	for(r=0;r<N_Z;r++){
		fscanf(fp, "%lf %lf", &temp, &(Zcrd[r]));
	}

	fclose(fp);

	return 0;
}

int Vec_Prep_Data::LoadReceivers(char *fname, int& _n_pointres, double (*(&_pointres))[3])
{
	FILE *fp=NULL;
	int i;

	if (_pointres) {delete [] _pointres; _pointres=NULL;}
	_n_pointres = 0;

	if((fp=fopen(fname, "r"))==0)
		Cannot_open_file(fname, "Vec_Prep_Data::LoadReceivers");

	fscanf(fp, "%ld", &_n_pointres);

	_pointres = new double[_n_pointres][3];
	if(_pointres == 0)
		Memory_allocation_error("_pointres", "Vec_Prep_Data::LoadReceivers");

	for(i=0; i<_n_pointres; i++)
		fscanf(fp, "%lf %lf %lf", &_pointres[i][0], &_pointres[i][1], &_pointres[i][2]);
	fclose(fp);

	return 0;
}
