#pragma once
#include "T_Brick.h"

class VecBlockSLAE
{
public:
	long *ig, *jg;

	double *di_block; // диагональные блоки
	double *gg_block; // внедиагональные блоки
	long *idi; // адреса хранения начал диагональных блоков
	long *ijg; // адреса хранения начал внедиагональных блоков

	double *pr; // глобальный вектор правой части

	long n_elem;    // число элементов в сетке
	long n_edges;   // число рёбер (всего)
	long n;         // размерность СЛАУ
	long nb;        // число блоков (блочная размерность СЛАУ)
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


