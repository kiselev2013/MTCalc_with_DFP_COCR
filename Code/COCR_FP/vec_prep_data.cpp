#include "stdafx.h" 
#include "vec_prep_data.h"

extern ofstream logfile;
//-----------------------------------------------------------------------------
Vec_Prep_Data::Vec_Prep_Data()
{
	nver = NULL;
	nvkat = NULL;
	xyz = NULL;
	mu3d = NULL;
	mu0 = NULL;
    sigma3d = NULL;
	sigma0 = NULL;
	l13d = NULL;
	dpr3d = NULL;
	dpr0 = NULL;
}
//-----------------------------------------------------------------------------
Vec_Prep_Data::~Vec_Prep_Data()
{
	if(mu0) {delete [] mu0; mu0=NULL;}
	if(xyz) {delete [] xyz; xyz=NULL;}
	if(nver) {delete [] nver; nver=NULL;}
	if(mu3d) {delete [] mu3d; mu3d=NULL;}
	if(l13d) {delete [] l13d; l13d=NULL;}
	if(nvkat) {delete [] nvkat; nvkat=NULL;}
	if(sigma0) {delete [] sigma0; sigma0=NULL;}
	if(sigma3d) {delete [] sigma3d; sigma3d=NULL;}
	if(dpr3d) {delete [] dpr3d; dpr3d=NULL;}
	if(dpr0) {delete [] dpr0; dpr0=NULL;}
}
//-----------------------------------------------------------------------------
// чтение сетки для МТЗ на векторных
//-----------------------------------------------------------------------------
int Vec_Prep_Data::Read_prep_data()
{
	In_Out R;
	FILE *fp=NULL;
	double temp1, temp2;
	int i, j;
	int max_material;
	int tmp;

	ifstream inf;

	// n_materials (максимальный номер материала в файле sig3d)
	if((fp=fopen("sig3d", "r"))==0)
		Cannot_open_file("sig3d", "Vec_Prep_Data::Read_prep_data");

	max_material = 0;
	while(!feof(fp))
	{
		fscanf(fp, "%d %lf %lf", &tmp, &temp1, &temp2);
		if(tmp > max_material)
			max_material = tmp;
	}
	n_materials = max_material;
	fclose(fp);

	// mu3d (mu3d, mu0)
	mu3d = new double[n_materials];
	if(mu3d == 0)
		Memory_allocation_error("mu3d", "Vec_Prep_Data::Read_prep_data");

	mu0 = new double[n_materials];
	if(mu0 == 0)
		Memory_allocation_error("mu0", "Vec_Prep_Data::Read_prep_data");

	if((fp=fopen("mu3d", "r"))==0)
		Cannot_open_file("mu3d", "Vec_Prep_Data::Read_prep_data");

	//while(!feof(fp))
	for (i=0; i<n_materials; i++)
	{
		int n_of_current_material;
		fscanf(fp, "%d", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf %lf",
			&mu3d[n_of_current_material], &mu0[n_of_current_material]);
		mu3d[n_of_current_material]*=MU_0;
		mu0[n_of_current_material]*=MU_0;
	}
	fclose(fp);

	// sig3d (sigma3d, sigma0)
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
		int n_of_current_material;
		fscanf(fp, "%d", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf %lf",
			&sigma3d[n_of_current_material], &sigma0[n_of_current_material]);
	}
	fclose(fp);

	// dpr3d (dpr3d, dpr0)
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
		int n_of_current_material;
		fscanf(fp, "%d", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf %lf",
			&dpr3d[n_of_current_material], &dpr0[n_of_current_material]);
		dpr3d[n_of_current_material]*=DPR_0;
		dpr0[n_of_current_material]*=DPR_0;
	}
	fclose(fp);

	// inftry.dat (kuzlov, kpar, kt1)
	R.Read_inftry("inftry.dat", &kuzlov, &kpar, &kt1);

	// nver.dat 
	nver = new int[kpar][14];
	if(nver == 0)
		Memory_allocation_error("nver", "Vec_Prep_Data::Read_prep_data");

	R.Read_Bin_File_Of_Long("nver.dat", (int*)nver, kpar, 14);
	for(i=0; i<kpar; i++)
		for(j=0; j<14; j++)
			nver[i][j]--;

	// nvkat.dat
	nvkat = new int[kpar];
	if(nvkat == 0)
		Memory_allocation_error("nvkat", "Vec_Prep_Data::Read_prep_data");

	R.Read_Bin_File_Of_Long("nvkat.dat", nvkat, kpar, 1);
	for(i=0; i<kpar; i++)
		nvkat[i]--;

	// xyz.dat
	xyz = new double[kuzlov][3];
	if(xyz == 0)
		Memory_allocation_error("xyz", "Vec_Prep_Data::Read_prep_data");

	R.Read_Bin_File_Of_Double("xyz.dat", (double*)xyz, kuzlov, 3);

	//l13d.dat
	l13d = new int[kt1];
	if(l13d == 0)
		Memory_allocation_error("l13d", "Vec_Prep_Data::Read_prep_data");

	R.Read_Bin_File_Of_Long("l13d.dat", l13d, kt1, 1);
	for(i=0; i<kt1; i++)
		l13d[i]--;

	return 0;
}
//------------------------------------------------------------------------
