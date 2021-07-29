// HarmVect3d.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "vec_prep_data.h"
#include "T_Mapping.h"
#include "ControlOMP.h"
#include "pcocr.h"
#include <chrono>

ofstream logfile;
ControlOMP omp;

const int size_i=sizeof(int);
const int size_d=sizeof(double);

int inputvec(int p_kpar,int &p_n,int &p_nc,int &p_vc,int (**ed)[25],int (**edges)[2],int **ig_t,int **jg_t,double **gg_t)
{
	int i,j;
	FILE *fp;
	ifstream inf;

	int mv;

	inf.open("tsize3d_.dat");
	if(!inf){
		printf("Error open file tsize3d_.dat");
		return 1;
	}
	inf>>p_vc;
	inf>>p_n;
	inf.close();
	inf.clear();

	p_nc=p_n-p_vc;
	if(!((*ig_t)=new int[p_n+1]))return 1;

	if (p_vc > 0)
	{
		if (!(fp = fopen("ig3d_.dat", "rb"))) {
			printf("Error open file ig3d.dat");
			return 1;
		}
		for (i = 0; i < p_nc; i++)(*ig_t)[i] = 0;
		fread(*ig_t + p_nc, size_i, p_vc + 1, fp);
		for (i = p_nc; i <= p_n; i++)(*ig_t)[i]--;
		fclose(fp);

		mv = (*ig_t)[p_n];

		if (!((*jg_t) = new int[mv]))return 1;
		if (!(fp = fopen("jg3d_.dat", "rb"))) {
			printf("Error open file jg3d.dat");
			return 1;
		}
		fread(*jg_t, size_i, mv, fp);
		for (i = 0; i < mv; i++)(*jg_t)[i]--;
		fclose(fp);

		if (!((*gg_t) = new double[mv]))return 1;
		if (!(fp = fopen("gg3d_.dat", "rb"))) {
			printf("Error open file gg3d.dat");
			return 1;
		}
		fread(*gg_t, size_d, mv, fp);
		fclose(fp);
	}

	if(!(fp=fopen("nodesforedges.dat","rb"))){
		printf("Error open file nodesforedges.dat");
		return 1;
	}

	if(!((*edges)=new int[p_n][2]))return 1;
	for(i=0;i<p_n;i++){
		fread((*edges)[i],size_i,2,fp);
		(*edges)[i][0]--;(*edges)[i][1]--;
	}
	fclose(fp);

	if(!(fp=fopen("edges.dat","rb"))){
		printf("Error open file edges.dat");
		return 1;
	}

	if(!((*ed)=new int[p_kpar][25]))return 1;
	for(i=0;i<p_kpar;i++){
		fread((*ed)[i],size_i,25,fp);
		for(j=0;j<25;j++)(*ed)[i][j]--;
	}
	fclose(fp);
	return 0;
}


//------------------------------------------------------------------------
int main()
{
	int i,e,j2;
	ifstream inf;
	ofstream ofp;
	int nproc;

	auto tProgramStart = std::chrono::high_resolution_clock::now();

	logfile.open("logharm3dCalc");

	omp.InitNumThreads();

	Vec_Prep_Data *d=NULL;
	if((d = new Vec_Prep_Data())==0) 
		throw logic_error("no memory for 3d mesh");

	// 0 - prosto raschet
	// 1 - s vikladkoi summarnogo polya po uzlam
	// 2 - dlya rascheta polya vliyaniya

	nproc=1;

	auto tNow = std::chrono::high_resolution_clock::now();
	logfile << std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(tNow - tProgramStart).count() << " : Reading TMatrix" << endl;
	inf.open("nthreads.txt");
	if(inf)
	{
		inf>>nproc;
		inf.close();
	}
	inf.clear();

	logfile << "nproc=" << nproc << endl;
	omp.SetNumThreads(nproc);
	logfile << "omp=" << omp.GetNumberOfThreads() << endl;

	T_Mapping_Vec *tmap=NULL;
	
	d->Read_prep_data();

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

	tNow = std::chrono::high_resolution_clock::now();
	logfile << std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(tNow - tProgramStart).count() << " : Reading SLAE" << endl;

	In_Out r;
	double *u=NULL;
	int n,n_edges_c,ig_n_1;
	double eps;
	int maxiter;

	inf.open("kuslau");
	if(!inf)
	{
		logfile<<"Error in open file "<<"kuslau"<<endl;
		cout<<"Error in open file "<<"kuslau"<<endl;
		exit(1);
	}
	inf>>n;
	inf>>eps;
	inf>>maxiter;
	inf.close();
	inf.clear();

	n_edges_c=tmap->n_c;

	int *ig=new int[n_edges_c+1];
	int *idi=new int[n_edges_c+1];
	double *pr=new double[n];

	r.Read_Bin_File_Of_Double("pr", pr, n, 1);

	r.Read_Bin_File_Of_Long("ig", ig, n_edges_c+1, 1);
	for(i=0; i<(n_edges_c+1); i++)
		ig[i]--;

	r.Read_Bin_File_Of_Long("idi", idi, n_edges_c+1, 1);
	for(i=0; i<(n_edges_c+1); i++)
		idi[i]--;

	ig_n_1=ig[n_edges_c];

	int *jg=new int[ig_n_1];
	int *ijg=new int[ig_n_1+1];
	double *di_block=new double[idi[n_edges_c]];


	r.Read_Bin_File_Of_Double("di", di_block, idi[n_edges_c], 1);

	
	r.Read_Bin_File_Of_Long("jg", jg, ig_n_1, 1);
	for(i=0; i<ig_n_1; i++)
		jg[i]--;

	r.Read_Bin_File_Of_Long("ijg", ijg, ig_n_1+1, 1);
	for(i=0; i<(ig_n_1+1); i++)
		ijg[i]--;

	double *gg_block=new double[ijg[ig_n_1]];
	r.Read_Bin_File_Of_Double("gg", gg_block, ijg[ig_n_1], 1);

	// решение СЛАУ

	if((u = new double[2*(tmap->n)])==0) Memory_allocation_error("u", "Loop_Harm_Vector_FEM");

	for (int i=0; i<n; i++)u[i]=0;
	int tsize3d;
	{
		double *sig = NULL;
		double *y_omp = NULL;
		int *ig_t=NULL;
		int *jg_t=NULL;
		double *gg_t=NULL;
		PCOCR pcocr;

		sig = new double[d->n_materials];
		y_omp = new double[n*omp.GetMaxThreads()];
		for (int i=0; i<d->n_materials; i++)
		{
			if (d->sigma3d[i] == 0)
				sig[i] = 0;
			else
				sig[i] = 1;
		}
		int n_nodes_c;
		if (r.Read_Long_From_Txt_File("tsize3d.dat", &tsize3d, false) != 0)
			tsize3d = 0;
		n_nodes_c = d->kuzlov - tsize3d;

		ig_t = new int[d->kuzlov + 1];

		if(tsize3d)
		{
			r.Read_Bin_File_Of_Long("ig3d.dat", &ig_t[n_nodes_c], tsize3d+1, 1);
			for(int i=0; i<n_nodes_c; i++) // дописываем единички в начало
			{
				ig_t[i] = 1;
			}
			for(int i=0; i<d->kuzlov+1; i++) // В Си нумерация с нуля
			{
				ig_t[i]--;
			}
		}
		else
		{
			for(int i=0; i<d->kuzlov+1; i++) // дописываем единички в начало
			{
				ig_t[i] = 0;
			}
		}
		// ig_t_n_1
		int ig_t_n_1 = ig_t[d->kuzlov];

		// gg3d.dat (gg_t)
		gg_t = new double[ig_t_n_1];

		if(tsize3d)
		{
			r.Read_Bin_File_Of_Double("gg3d.dat", gg_t, ig_t_n_1, 1);
		}

		// jg3d.dat (jg_t)
		jg_t = new int[ig_t_n_1];

		if(tsize3d)
		{
			r.Read_Bin_File_Of_Long("jg3d.dat", jg_t, ig_t_n_1, 1);
		}
		for(int i=0; i<ig_t_n_1; i++)
		{
			jg_t[i]--;
		}

		int *is_node_bound=NULL;

		if((is_node_bound = new int[d->kuzlov])==0) 
			Memory_allocation_error("is_node_bound","Bound_cond_vec_harm::Make_list_of_bound_edges_harm");

		for(int i=0; i<d->kuzlov; i++)
		{
			is_node_bound[i] = 0;
		}

		for(int i=0; i<d->kt1; i++)
		{
			is_node_bound[d->l13d[i]] = 1;
		}
	

		tNow = std::chrono::high_resolution_clock::now();
		logfile << std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(tNow - tProgramStart).count() << " : Solving SLAE" << endl;

		pcocr.PCOCR_2x2_Folded(n, ig, jg, idi, ijg, di_block, gg_block, pr, u, eps, maxiter, y_omp, d->kpar, tmap->n_c, n_nodes_c, d->xyz, d->nvkat, d->nver, sig, tmap->edges, ig_t, jg_t, gg_t, is_node_bound);


		delete [] is_node_bound;
		delete [] sig;
		delete [] y_omp;
		delete [] ig_t;
		delete [] jg_t;
		delete [] gg_t;
	}

	tNow = std::chrono::high_resolution_clock::now();
	logfile << std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(tNow - tProgramStart).count() << " : Writing result" << endl;
	// запись решения в файл
	r.Write_Bin_File_Of_Double("v3.dat", u, n, 1);

	if (tsize3d > 0)
	{
		for (e = tmap->n_c; e < tmap->n; e++)
		{
			u[2 * e] = 0.0;
			u[2 * e + 1] = 0.0;
			for (j2 = tmap->ig_t[e]; j2 <= tmap->ig_t[e + 1] - 1; j2++)
			{
				u[2 * e] += u[2 * (tmap->jg_t[j2])] * tmap->gg_t[j2];
				u[2 * e + 1] += u[2 * (tmap->jg_t[j2]) + 1] * tmap->gg_t[j2];
			}
		}
		r.Write_Bin_File_Of_Double("v3.upk", u, tmap->n, 2);
	}

	if(d) {delete d; d=NULL;}
	if(u) {delete [] u; u=NULL;}
	if(tmap) {delete tmap; tmap=NULL;}


	tNow = std::chrono::high_resolution_clock::now();
	logfile << std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> >>(tNow - tProgramStart).count() << " : All done" << endl;

	logfile.close();
	logfile.clear();

	return 0;
}

