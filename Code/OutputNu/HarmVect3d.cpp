/**
 * GENERAL REMARKS
 *
 *  This code is freely available under the following conditions:
 *
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer.
 *
 *              MTCalc_with_DFP_COCR
 *  Output 3D EM field with smoothing
 *
 *  Written by Prof. Marina G. Persova and Ph.D. Petr A. Domnikov 
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  mpersova@mail.ru
 *  Version 1.5 December 17, 2020
*/

#include "stdafx.h"
#include "AbstractFEM.h"
#include "T_Mapping.h"
#include "give_out_vec_mt.h"
#include "vec_prep_data.h"
#include "ListForLine.h"
#include "in_out.h"

ofstream logfile;

const int size_i=sizeof(int);
const int size_d=sizeof(double);

//-----------------------------------------------------------
// Reading an information about edges
//-----------------------------------------------------------
int inputvec(long p_kpar,long &p_n,long &p_nc,long &p_vc,long (**ed)[25],long (**edges)[2],long **ig_t,long **jg_t,double **gg_t)
{
	int i,j;
	FILE *fp;
	ifstream inf;
	int mv;

	inf.open("tsize3d_.dat");
	if(!inf){
		printf("Error open file tsize3d.dat");
		return 1;
	}
	inf>>p_vc;
	inf>>p_n;
	inf.close();
	inf.clear();

	p_nc=p_n-p_vc;
	if(!((*ig_t)=new long[p_n+1]))return 1;

	if(!(fp=fopen("ig3d_.dat","rb"))){
		printf("Error open file ig3d.dat");
		return 1;
	}
	for(i=0;i<p_nc;i++)(*ig_t)[i]=0;
	fread(*ig_t+p_nc,size_i,p_vc+1,fp);
	for(i=p_nc;i<=p_n;i++)(*ig_t)[i]--;
	fclose(fp);

	mv=(*ig_t)[p_n];

	if(!((*jg_t)=new long[mv]))return 1;
	if(!(fp=fopen("jg3d_.dat","rb"))){
		printf("Error open file jg3d.dat");
		return 1;
	}
	fread(*jg_t,size_i,mv,fp);
	for(i=0;i<mv;i++)(*jg_t)[i]--;
	fclose(fp);

	if(!((*gg_t)=new double[mv]))return 1;
	if(!(fp=fopen("gg3d_.dat","rb"))){
		printf("Error open file gg3d.dat");
		return 1;
	}
	fread(*gg_t,size_d,mv,fp);
	fclose(fp);


	if(!(fp=fopen("nodesforedges.dat","rb"))){
		printf("Error open file nodesforedges.dat");
		return 1;
	}

	if(!((*edges)=new long[p_n][2]))return 1;
	for(i=0;i<p_n;i++){
		fread((*edges)[i],size_i,2,fp);
		(*edges)[i][0]--;(*edges)[i][1]--;
	}
	fclose(fp);

	if(!(fp=fopen("edges.dat","rb"))){
		printf("Error open file edges.dat");
		return 1;
	}

	if(!((*ed)=new long[p_kpar][25]))return 1;
	for(i=0;i<p_kpar;i++){
		fread((*ed)[i],size_i,25,fp);
		for(j=0;j<25;j++)(*ed)[i][j]--;
	}
	fclose(fp);
	return 0;
}
//-----------------------------------------------------------
// Main function of module
//-----------------------------------------------------------
int main()
{
	int i;
	ifstream inf;
	double tmpd;

	logfile.open("logharm3dOutput");

	Vec_Prep_Data *d=NULL;
	if((d = new Vec_Prep_Data())==0) 
		throw logic_error("no memory for 3d mesh");

	d->tasktype=0;

	T_Mapping_Vec *tmap=NULL;
	In_Out r;
	
	d->Read_prep_data();

	double (*xyz0)[3];
	xyz0=NULL;
	d->xyzt=NULL;
	
	inf.open("xyz0.dat",ios::binary);
	if(inf)
	{
		xyz0=new double[d->kuzlov][3];
		for(i=0;i<d->kuzlov;i++)
		{
			inf>xyz0[i][0]>xyz0[i][1]>xyz0[i][2];
		}
		inf.close();
	}
	inf.clear();

	if(!d->fdirect)
	{
		d->Read_mtz_1d();
	}
	else
	{
		inf.open("alfa");
		if(!inf)
		{
			logfile<<"Error in open file "<<"alfa"<<endl;
			logfile.close();
			return 1;
		}
		inf>>tmpd;
		inf.close();
		inf.clear();
		if (fabs(1.0 - tmpd)<0.01)
		{
			d->alfa = 1; 
		} 
		else
		{
			d->alfa = 0;
		}
	}

	if((tmap = new T_Mapping_Vec(d->nver, d->xyz, d->kuzlov, d->kpar))==0)
	{
		if(d) {delete d; d=NULL;}
		throw logic_error("no memory for tmap");
	}
	
	if(inputvec(tmap->kpar,tmap->n,tmap->n_c,tmap->n_dc,&(tmap->ed),&(tmap->edges),&(tmap->ig_t),&(tmap->jg_t),&(tmap->gg_t)))
	{
		logfile<<"Error in inputvecmain"<<endl;
		logfile.close();
		return 1;
	}

	double *u=NULL;

	if((u = new double[2*(tmap->n)])==0) Memory_allocation_error("u", "main");
	
	r.Read_Bin_File_Of_Double("v3.dat", u, tmap->n_c, 2);

	Give_out_vec_mt *give_out=NULL;
	
	int e,j2;
	for(e=tmap->n_c;e<tmap->n;e++)
	{
		u[2*e] = 0.0;
		u[2*e+1] = 0.0;
		for(j2=tmap->ig_t[e]; j2<=tmap->ig_t[e+1]-1; j2++)
		{
			u[2*e] += u[2*(tmap->jg_t[j2])]*tmap->gg_t[j2];
			u[2*e+1] += u[2*(tmap->jg_t[j2])+1]*tmap->gg_t[j2];
		}
	}

	give_out = new Give_out_vec_mt(d, tmap, u);
	d->xyzt=d->xyz;
	if(xyz0)
	{
		d->xyz=xyz0;
	}
	give_out->resultantB = new OutputResultant3d(give_out,vtWithoutDiscontinuity);
	give_out->resultantB->Prepare();
	give_out->resultantA = new OutputResultant3d(give_out,vtWithDiscontinuity);
	give_out->resultantA->Prepare();
	if(xyz0)
	{
		d->xyz=d->xyzt;
		d->xyzt=NULL;
	}
	give_out->Give_out_on_hex(true);
	give_out->resultantB->StopSolvers();
	delete give_out->resultantB;
	give_out->resultantB = NULL;
	give_out->resultantA->StopSolvers();
	delete give_out->resultantA;
	give_out->resultantA = NULL;
	delete give_out;

	if(xyz0){delete [] xyz0;xyz0=NULL;}

	if(u) {delete [] u; u=NULL;}
	if(d) {delete d; d=NULL;}
	if(tmap) {delete tmap; tmap=NULL;}

	logfile.close();

	return 0;
}
