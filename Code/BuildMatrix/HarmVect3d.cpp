#include "stdafx.h"
#include "bound_cond_vec_harm.h"
#include "t_global_slae.h"
#include "T_Portrait.h"

ofstream logfile;

const int size_i=sizeof(int);
const int size_d=sizeof(double);

int inputvec(long p_kpar,long &p_n,long (**ed)[25],long (**edges)[2])
{
	int i,j,p_vc;
	FILE *fp;
	ifstream inf;

	inf.open("tsize3d_.dat");
	if(!inf){
		printf("Error open file tsize3d_.dat");
		return 1;
	}
	inf>>p_vc;
	inf>>p_n;
	inf.close();
	inf.clear();

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

int main()
{
	ifstream inf;
	ofstream ofp;

	long n;				// всего базисных функций
	long (*edges)[2];	// рЄбра, заданные 2-м€ вершинами
	long (*ed)[25];		// €чейки перечисленные своими рЄбрами + доп. информаци€
	
	logfile.open("logharm3d");

	logfile<<"Start program "<<(int)clock()<<'\n';

	Vec_Prep_Data *d=NULL;
	if((d = new Vec_Prep_Data())==0) 
		throw logic_error("no memory for 3d mesh");

	d->Read_prep_data();
	d->Read_mtz_1d();

	if(inputvec(d->kpar,n,&(ed),&(edges)))
	{
		logfile<<"Error in inputvecmain"<<endl;
		logfile.close();
		return 1;
	}

	T_Portrait *p=NULL;
	VecBlockSLAE *slae=NULL;
	Bound_cond_vec_harm *bc=NULL;

	if((p = new T_Portrait((long*)ed,n,d->kpar))==0)
		Memory_allocation_error("p", "Loop_Harm_Vector_FEM");

	logfile<<"Begin build portrait "<<(int)clock()<<'\n';

	p->Gen_Portrait();
	p->Gen_idi_ijg(d->nvkat, d->nver);

	logfile<<"Finish build portrait "<<(int)clock()<<'\n';
	logfile<<"Start assembling SLAE "<<(int)clock()<<'\n';

	if((slae = new VecBlockSLAE(p->ig, p->jg, p->idi, p->ijg, d->kpar, n, d->xyz, d->nver, ed, edges, 
		d->nvkat, d->nu,  d->mu3d, d->mu0, d->sigma3d, d->sigma0, d->dpr3d, d->dpr0,
		d->alfa, d->n_1d, d->z_1d, d->usin, d->ucos))==0)
		Memory_allocation_error("slae", "Loop_Harm_Vector_FEM");

	slae->AsmBlockSLAE(d);

	if((bc = new Bound_cond_vec_harm(d->kuzlov, n, d->kt1, d->l13d, edges, d->kpar))==0)
		Memory_allocation_error("bc", "Loop_Harm_Vector_FEM");

	bc->SetHomogenDirichletCond(p->ig, p->jg, p->idi, p->ijg, slae->di_block, slae->gg_block, slae->pr);

	logfile<<"Finish assembling SLAE "<<(int)clock()<<'\n';

	slae->WriteBlockSLAEtoFiles();

	logfile<<"Finish program "<<(int)clock()<<'\n';

	logfile.close();
	logfile.clear();

	return 0;
}
