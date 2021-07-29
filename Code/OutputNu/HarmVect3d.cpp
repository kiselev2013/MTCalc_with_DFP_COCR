#include "stdafx.h"
#include "AbstractFEM.h"
#include "T_Mapping.h"
#include "give_out_vec_mt.h"
#include "vec_prep_data.h"
#include "GeoPrepDocSettings.h"
#include "ListForLine.h"
#include "in_out.h"
#include "ControlOMP.h"

GeoPrepDocSettings GPDocSettings;

ofstream logfile;

ControlOMP omp;

const int size_i=sizeof(int);
const int size_d=sizeof(double);

struct Point3D
{
	double x,y,z;
};

extern void GetNodeSolutionByLocalcCoords(double *lc,double *q,double &Val);

void GetGlobalCoordinates(Point3D *HexPnt,double *lc,double *gc)
{
	double ksi,eta,phi;

	ksi=lc[0];
	eta=lc[1];
	phi=lc[2];

	gc[0]=	HexPnt[0].x*(1-ksi)*(1-eta)*(1-phi)+HexPnt[1].x*(ksi)*(1-eta)*(1-phi)+
			HexPnt[2].x*(1-ksi)*(eta)*(1-phi)+HexPnt[3].x*(ksi)*(eta)*(1-phi)+
			HexPnt[4].x*(1-ksi)*(1-eta)*(phi)+HexPnt[5].x*(ksi)*(1-eta)*(phi)+
			HexPnt[6].x*(1-ksi)*(eta)*(phi)+HexPnt[7].x*(ksi)*(eta)*(phi);
	gc[1]=	HexPnt[0].y*(1-ksi)*(1-eta)*(1-phi)+HexPnt[1].y*(ksi)*(1-eta)*(1-phi)+
			HexPnt[2].y*(1-ksi)*(eta)*(1-phi)+HexPnt[3].y*(ksi)*(eta)*(1-phi)+
			HexPnt[4].y*(1-ksi)*(1-eta)*(phi)+HexPnt[5].y*(ksi)*(1-eta)*(phi)+
			HexPnt[6].y*(1-ksi)*(eta)*(phi)+HexPnt[7].y*(ksi)*(eta)*(phi);
	gc[2]=	HexPnt[0].z*(1-ksi)*(1-eta)*(1-phi)+HexPnt[1].z*(ksi)*(1-eta)*(1-phi)+
			HexPnt[2].z*(1-ksi)*(eta)*(1-phi)+HexPnt[3].z*(ksi)*(eta)*(1-phi)+
			HexPnt[4].z*(1-ksi)*(1-eta)*(phi)+HexPnt[5].z*(ksi)*(1-eta)*(phi)+
			HexPnt[6].z*(1-ksi)*(eta)*(phi)+HexPnt[7].z*(ksi)*(eta)*(phi);

}

void FindLocalCoordinates(Point3D &R,Point3D *HexPnt,double *lc)
{
	double EpsForFindLocalCoord = 1e-2;
	double CoeffForDiv = 1.2;
	int MaxDeep = 10;


	int p,t,m,deep,crd[3][2]={{1,2},{0,2},{0,1}};
	double LocalCoord[2][3],CentGlob[3],CentLoc[3],DopLoc[3];
	double dist,disc,h[3];
	const double ods=0.16666666666666666;
	const double hds=0.33333333333333333;
	deep=0;
	LocalCoord[0][0]=0.0;LocalCoord[0][1]=0.0;LocalCoord[0][2]=0.0;
	LocalCoord[1][0]=1.0;LocalCoord[1][1]=1.0;LocalCoord[1][2]=1.0;
	CentLoc[0]=0.5;CentLoc[1]=0.5;CentLoc[2]=0.5;
	GetGlobalCoordinates(HexPnt,CentLoc,CentGlob);
	disc=sqrt(
		(R.x-CentGlob[0])*(R.x-CentGlob[0])+
		(R.y-CentGlob[1])*(R.y-CentGlob[1])+
		(R.z-CentGlob[2])*(R.z-CentGlob[2])
		);
	h[0]=LocalCoord[1][0]-LocalCoord[0][0];
	h[1]=LocalCoord[1][1]-LocalCoord[0][1];
	h[2]=LocalCoord[1][2]-LocalCoord[0][2];
	do{
		if(disc<EpsForFindLocalCoord)break;

		for(m=0;m<3;m++){

			CentLoc[m]=LocalCoord[0][m]+ods*h[m];
			
			GetGlobalCoordinates(HexPnt,CentLoc,CentGlob);
			dist=sqrt(
					(R.x-CentGlob[0])*(R.x-CentGlob[0])+
					(R.y-CentGlob[1])*(R.y-CentGlob[1])+
					(R.z-CentGlob[2])*(R.z-CentGlob[2])
					);

			if(disc<dist){
				CentLoc[m]=LocalCoord[1][m]-ods*h[m];
			
				GetGlobalCoordinates(HexPnt,CentLoc,CentGlob);
				dist=sqrt(
						(R.x-CentGlob[0])*(R.x-CentGlob[0])+
						(R.y-CentGlob[1])*(R.y-CentGlob[1])+
						(R.z-CentGlob[2])*(R.z-CentGlob[2])
						);
				if(dist<disc){
					disc=dist;
					LocalCoord[0][m]=LocalCoord[1][m]-hds*h[m];
				}
				else{
					LocalCoord[0][m]=LocalCoord[0][m]+hds*h[m];
					LocalCoord[1][m]=LocalCoord[1][m]-hds*h[m];
				}
			}
			else{
				disc=dist;
				LocalCoord[1][m]=LocalCoord[0][m]+hds*h[m];				
			}
			h[m]=LocalCoord[1][m]-LocalCoord[0][m];
			CentLoc[m]=0.5*(LocalCoord[0][m]+LocalCoord[1][m]);
		}
		deep++;
	}while(deep<MaxDeep);
	
	// Допоиск
	if(deep==MaxDeep){
		DopLoc[0]=CentLoc[0];
		DopLoc[1]=CentLoc[1];
		DopLoc[2]=CentLoc[2];
		do{
			t=0;
			for(m=0;m<3;m++){
				do{
					p=0;
					DopLoc[m]=CentLoc[m]-h[m];				
					if(DopLoc[m]>0){
						GetGlobalCoordinates(HexPnt,DopLoc,CentGlob);
						dist=sqrt(
								(R.x-CentGlob[0])*(R.x-CentGlob[0])+
								(R.y-CentGlob[1])*(R.y-CentGlob[1])+
								(R.z-CentGlob[2])*(R.z-CentGlob[2])
								);
						if(dist<disc){
							disc=dist;
							CentLoc[m]=DopLoc[m];
							t=p=1;
						}
					}
				}while(p);
				do{
					p=0;
					DopLoc[m]=CentLoc[m]+h[m];				
					if(DopLoc[m]<1){
						GetGlobalCoordinates(HexPnt,DopLoc,CentGlob);
						dist=sqrt(
								(R.x-CentGlob[0])*(R.x-CentGlob[0])+
								(R.y-CentGlob[1])*(R.y-CentGlob[1])+
								(R.z-CentGlob[2])*(R.z-CentGlob[2])
								);
						if(dist<disc){
							disc=dist;
							CentLoc[m]=DopLoc[m];
							t=p=1;
						}
					}
				}while(p);
			}
		}while(t);
	}

	GetGlobalCoordinates(HexPnt,CentLoc,CentGlob);

	if(disc>EpsForFindLocalCoord){
		logfile<<"Reciver Eps Fail:"<<endl;
		logfile<<R.x<<'\t'<<R.y<<'\t'<<R.z<<endl;
		logfile<<CentGlob[0]<<'\t'<<CentGlob[1]<<'\t'<<CentGlob[2]<<endl;
		logfile<<sqrt(R.x*R.x+R.y*R.y+R.z*R.z)<<endl;
		logfile<<sqrt((R.x-CentGlob[0])*(R.x-CentGlob[0])+
				   (R.y-CentGlob[1])*(R.y-CentGlob[1])+
				   (R.z-CentGlob[2])*(R.z-CentGlob[2]))<<endl;
		logfile<<endl;
	}
	
	lc[0]=CentLoc[0];
	lc[1]=CentLoc[1];
	lc[2]=CentLoc[2];
}

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
//------------------------------------------------------------------------
void OutEnormHarm(Vec_Prep_Data *d,T_Mapping_Vec *tmap,double *u)
{
	int nmat;
	vector<int> ElemsSizeForNodes,ElemsForNodes,ElemsIgNodes;
	vector<int> MatsForNodes,MatsIgNodes,MatsSizeForNodes;
	vector<int> ElemsForMats,ElemsIgMats,ElemsSizeForMats;
	vector<double> RegularLocalCoord;

	int i,j,k,n,m,nigsz,li,lj,lk,ln,lp,lt,kuzlov_reg,kpar_reg,nx,ny,nz,kx,ky,kz,e;
	double lc[3];

	nx=d->N_X;
	ny=d->N_Y;
	nz=d->N_Z;
	kx=nx-1;
	ky=ny-1;
	kz=nz-1;
	kuzlov_reg=nx*ny*nz;
	kpar_reg=kx*ky*kz;

	ListOfNodesForLine nodesLine;

	ElemsSizeForNodes.resize(kuzlov_reg);
	ElemsIgNodes.resize(kuzlov_reg+1);
	for(i=0;i<kuzlov_reg;i++)ElemsSizeForNodes[i]=0;

	for(lk=0;lk<kz;lk++)
	{
		for(lj=0;lj<ky;lj++)
		{
			for(li=0;li<kx;li++)
			{
				// Одни и теже нерегулярные элементы посчитаются несколько раз
				ElemsSizeForNodes[nx*(ny*lk+lj)+li]++;
				ElemsSizeForNodes[nx*(ny*lk+lj)+(li+1)]++;
				ElemsSizeForNodes[nx*(ny*lk+(lj+1))+li]++;
				ElemsSizeForNodes[nx*(ny*lk+(lj+1))+(li+1)]++;
				ElemsSizeForNodes[nx*(ny*(lk+1)+lj)+li]++;
				ElemsSizeForNodes[nx*(ny*(lk+1)+lj)+(li+1)]++;
				ElemsSizeForNodes[nx*(ny*(lk+1)+(lj+1))+li]++;
				ElemsSizeForNodes[nx*(ny*(lk+1)+(lj+1))+(li+1)]++;
			}
		}
	}

	nigsz=0;
	for(i=0;i<kuzlov_reg;i++)nigsz+=ElemsSizeForNodes[i];
	ElemsForNodes.resize(nigsz);

	ElemsIgNodes[i]=0;
	for(i=0;i<kuzlov_reg;i++)
	{
		ElemsIgNodes[i+1]=ElemsIgNodes[i]+ElemsSizeForNodes[i];
		ElemsSizeForNodes[i]=0;
	}

	for(lk=0;lk<kz;lk++)
	{
		for(lj=0;lj<ky;lj++)
		{
			for(li=0;li<kx;li++)
			{
				i=d->regular[kx*(ky*lk+lj)+li];

				nigsz=nx*(ny*lk+lj)+li;
				ElemsForNodes[ElemsIgNodes[nigsz]+ElemsSizeForNodes[nigsz]]=i;
				ElemsSizeForNodes[nigsz]++;

				nigsz=nx*(ny*lk+lj)+(li+1);
				ElemsForNodes[ElemsIgNodes[nigsz]+ElemsSizeForNodes[nigsz]]=i;
				ElemsSizeForNodes[nigsz]++;

				nigsz=nx*(ny*lk+(lj+1))+li;
				ElemsForNodes[ElemsIgNodes[nigsz]+ElemsSizeForNodes[nigsz]]=i;
				ElemsSizeForNodes[nigsz]++;

				nigsz=nx*(ny*lk+(lj+1))+(li+1);
				ElemsForNodes[ElemsIgNodes[nigsz]+ElemsSizeForNodes[nigsz]]=i;
				ElemsSizeForNodes[nigsz]++;

				nigsz=nx*(ny*(lk+1)+lj)+li;
				ElemsForNodes[ElemsIgNodes[nigsz]+ElemsSizeForNodes[nigsz]]=i;
				ElemsSizeForNodes[nigsz]++;

				nigsz=nx*(ny*(lk+1)+lj)+(li+1);
				ElemsForNodes[ElemsIgNodes[nigsz]+ElemsSizeForNodes[nigsz]]=i;
				ElemsSizeForNodes[nigsz]++;

				nigsz=nx*(ny*(lk+1)+(lj+1))+li;
				ElemsForNodes[ElemsIgNodes[nigsz]+ElemsSizeForNodes[nigsz]]=i;
				ElemsSizeForNodes[nigsz]++;

				nigsz=nx*(ny*(lk+1)+(lj+1))+(li+1);
				ElemsForNodes[ElemsIgNodes[nigsz]+ElemsSizeForNodes[nigsz]]=i;
				ElemsSizeForNodes[nigsz]++;
			}
		}
	}

	// Убираем повторы
	for(i=0;i<kuzlov_reg;i++)
	{
		// Заменяем повторы на -1
		for(lj=ElemsIgNodes[i];lj<ElemsIgNodes[i+1];lj++)
		{
			if(ElemsForNodes[lj]!=-1)
			{
				for(lk=lj+1;lk<ElemsIgNodes[i+1];lk++)
				{
					if(ElemsForNodes[lj]==ElemsForNodes[lk])ElemsForNodes[lk]=-1;
				}
			}
		}
		// Уплотняем и подсчитываем колличество
		lk=ElemsIgNodes[i]+1;
		for(lj=ElemsIgNodes[i]+1;lj<ElemsIgNodes[i+1];lj++)
		{
			if(ElemsForNodes[lj]!=-1)
			{
				ElemsForNodes[lk]=ElemsForNodes[lj];
				lk++;
			}
		}
		ElemsSizeForNodes[i]=lk-ElemsIgNodes[i];
		for(;lk<ElemsIgNodes[i+1];lk++)
		{
			ElemsForNodes[lk]=-1;
		}
	}

	// Сдвигаем объекты и исправляем ig
	nigsz=ElemsSizeForNodes[0];
	for(i=1;i<kuzlov_reg;i++)
	{
		lj=nigsz;
		for(j=0;j<ElemsSizeForNodes[i];j++)
		{
			ElemsForNodes[nigsz]=ElemsForNodes[ElemsIgNodes[i]+j];
			nigsz++;
		}
		ElemsIgNodes[i]=lj;
	}
	ElemsIgNodes[kuzlov_reg]=nigsz;
	ElemsForNodes.resize(nigsz);

	nodesLine.item.resize(d->kuzlov);

	// Ребра для узлов (с каждым узлом может быть связано 6 обычных ребер и 6 терминальных)
	vector<int> efnd[12];
	vector<int> efndsz;
	int cnmat,crmat[8],ckmat[8];
	double FldValSin[64][3],FldValCos[64][3];

	efndsz.resize(d->kuzlov);
	for(k=0;k<d->kuzlov;k++)efndsz[k]=0;

	for(i=0;i<12;i++)
	{
		efnd[i].resize(d->kuzlov);
		for(k=0;k<d->kuzlov;k++)efnd[i][k]=-1;
	}

	for(i=0;i<tmap->n;i++)
	{
		for(j=0;j<2;j++)
		{
			n=tmap->edges[i][j];
			efnd[efndsz[n]][n]=i;
			efndsz[n]++;
		}
	}

	for(i=0;i<d->kuzlov;i++)
	{
		cnmat=0;
		for(k=0;k<8;k++)ckmat[k]=0;
		for(j=0;j<efndsz[i];j++)
		{
			n=efnd[j][i];
			m=(int)d->EnForLine[n].size();
			for(k=0;k<m;k++)
			{
				for(lp=cnmat-1;lp>=0;lp--)
				{
					if(crmat[lp]==d->EnForLine[n][k].mtr)break;
				}
				if(lp==-1)
				{
					crmat[cnmat]=d->EnForLine[n][k].mtr;
					lk=8*cnmat;
					ckmat[cnmat]=1;
					cnmat++;
				}
				else
				{
					lk=8*lp+ckmat[lp];
					ckmat[lp]++;
				}
				FldValSin[lk][0]=d->EnForLine[n][k].es[0];
				FldValSin[lk][1]=d->EnForLine[n][k].es[1];
				FldValSin[lk][2]=d->EnForLine[n][k].es[2];
				FldValCos[lk][0]=d->EnForLine[n][k].ec[0];
				FldValCos[lk][1]=d->EnForLine[n][k].ec[1];
				FldValCos[lk][2]=d->EnForLine[n][k].ec[2];
			}
		}
		nodesLine.item[i].resize(cnmat);
		for(lj=0;lj<cnmat;lj++)
		{
			nodesLine.item[i][lj].material=crmat[lj];
			nodesLine.item[i][lj].val[0][0]=0.0;
			nodesLine.item[i][lj].val[1][0]=0.0;
			nodesLine.item[i][lj].val[2][0]=0.0;
			nodesLine.item[i][lj].val[0][1]=0.0;
			nodesLine.item[i][lj].val[1][1]=0.0;
			nodesLine.item[i][lj].val[2][1]=0.0;
			for(li=0;li<ckmat[lj];li++)
			{
				lk=8*lj+li;
				nodesLine.item[i][lj].val[0][0]+=FldValSin[lk][0];
				nodesLine.item[i][lj].val[1][0]+=FldValSin[lk][1];
				nodesLine.item[i][lj].val[2][0]+=FldValSin[lk][2];
				nodesLine.item[i][lj].val[0][1]+=FldValCos[lk][0];
				nodesLine.item[i][lj].val[1][1]+=FldValCos[lk][1];
				nodesLine.item[i][lj].val[2][1]+=FldValCos[lk][2];
			}
			nodesLine.item[i][lj].val[0][0]/=ckmat[lj];
			nodesLine.item[i][lj].val[1][0]/=ckmat[lj];
			nodesLine.item[i][lj].val[2][0]/=ckmat[lj];
			nodesLine.item[i][lj].val[0][1]/=ckmat[lj];
			nodesLine.item[i][lj].val[1][1]/=ckmat[lj];
			nodesLine.item[i][lj].val[2][1]/=ckmat[lj];
		}
	}

	logfile<<"Calculate local coords for regular mesh"<<endl;
	Point3D R,CntHex[8];
	RegularLocalCoord.resize(3*nigsz);
	for(lk=0;lk<nz;lk++)
	{
		for(lj=0;lj<ny;lj++)
		{
			for(li=0;li<nx;li++)
			{
				i=nx*(ny*lk+lj)+li;
				for(lp=0;lp<ElemsSizeForNodes[i];lp++)
				{
					j=ElemsForNodes[ElemsIgNodes[i]+lp];
					if(!(d->nver[j][13]>30))
					{
						RegularLocalCoord[3*(ElemsIgNodes[i]+lp)]=(d->Xcrd[li]-d->xyz[d->nver[j][0]][0])/(d->xyz[d->nver[j][7]][0]-d->xyz[d->nver[j][0]][0]);
						RegularLocalCoord[3*(ElemsIgNodes[i]+lp)+1]=(d->Ycrd[lj]-d->xyz[d->nver[j][0]][1])/(d->xyz[d->nver[j][7]][1]-d->xyz[d->nver[j][0]][1]);
						RegularLocalCoord[3*(ElemsIgNodes[i]+lp)+2]=(d->Zcrd[lk]-d->xyz[d->nver[j][0]][2])/(d->xyz[d->nver[j][7]][2]-d->xyz[d->nver[j][0]][2]);
					}
					else
					{
						for(ln=0; ln<8; ln++)
						{
							CntHex[ln].x=d->xyz[d->nver[j][ln]][0];
							CntHex[ln].y=d->xyz[d->nver[j][ln]][1];
							CntHex[ln].z=d->xyz[d->nver[j][ln]][2];
						}
						R.x=d->Xcrd[li];
						R.y=d->Ycrd[lj];
						R.z=d->Zcrd[lk];

						FindLocalCoordinates(R,CntHex,lc);

						RegularLocalCoord[3*(ElemsIgNodes[i]+lp)]=lc[0];
						RegularLocalCoord[3*(ElemsIgNodes[i]+lp)+1]=lc[1];
						RegularLocalCoord[3*(ElemsIgNodes[i]+lp)+2]=lc[2];
					}
				}
			}
		}
	}

	logfile<<"Nodes to Materials"<<endl;
	// Материалы для узлов
	MatsSizeForNodes.resize(kuzlov_reg);
	MatsIgNodes.resize(kuzlov_reg+1);
	MatsForNodes.resize(nigsz);

	// Заменяем элементы их материалами
	for(i=0;i<kuzlov_reg;i++)
	{
		MatsSizeForNodes[i]=ElemsSizeForNodes[i];
		MatsIgNodes[i]=ElemsIgNodes[i];
	}
	MatsIgNodes[kuzlov_reg]=ElemsIgNodes[kuzlov_reg];
	for(i=0;i<nigsz;i++)MatsForNodes[i]=d->nvkat[ElemsForNodes[i]];

	// Убираем повторы
	for(i=0;i<kuzlov_reg;i++)
	{
		// Заменяем повторы на -1
		for(lj=MatsIgNodes[i];lj<MatsIgNodes[i+1];lj++)
		{
			if(MatsForNodes[lj]!=-1)
			{
				for(lk=lj+1;lk<MatsIgNodes[i+1];lk++)
				{
					if(MatsForNodes[lj]==MatsForNodes[lk])MatsForNodes[lk]=-1;
				}
			}
		}
		// Уплотняем и подсчитываем колличество
		lk=MatsIgNodes[i]+1;
		for(lj=MatsIgNodes[i]+1;lj<MatsIgNodes[i+1];lj++)
		{
			if(MatsForNodes[lj]!=-1)
			{
				MatsForNodes[lk]=MatsForNodes[lj];
				lk++;
			}
		}
		MatsSizeForNodes[i]=lk-MatsIgNodes[i];
		for(;lk<MatsIgNodes[i+1];lk++){
			MatsForNodes[lk]=-1;
		}
	} 

	// Сдвигаем объекты и исправляем ig
	nigsz=MatsSizeForNodes[0];
	for(i=1;i<kuzlov_reg;i++)
	{
		lj=nigsz;
		for(j=0;j<MatsSizeForNodes[i];j++)
		{
			MatsForNodes[nigsz]=MatsForNodes[MatsIgNodes[i]+j];
			nigsz++;
		}
		MatsIgNodes[i]=lj;
	}
	MatsIgNodes[kuzlov_reg]=nigsz;
	MatsForNodes.resize(nigsz);

	logfile<<"Materials to Elements"<<endl;
	// Формируем список элементов для материалов
	nmat=0;
	for(i=0;i<d->kpar;i++)
	{
		if(d->nvkat[i]+1>nmat)nmat=d->nvkat[i]+1;
	}
	
	ElemsSizeForMats.resize(nmat);
	ElemsIgMats.resize(nmat+1);
	
	for(i=0;i<nmat;i++)ElemsSizeForMats[i]=0;
	for(i=0;i<d->kpar;i++)ElemsSizeForMats[d->nvkat[i]]++;

	nigsz=0;
	for(i=0;i<nmat;i++)nigsz+=ElemsSizeForMats[i];
	ElemsForMats.resize(nigsz);

	ElemsIgMats[i]=0;
	for(i=0;i<nmat;i++)
	{
		ElemsIgMats[i+1]=ElemsIgMats[i]+ElemsSizeForMats[i];
		ElemsSizeForMats[i]=0;
	}

	for(i=0;i<d->kpar;i++)
	{
		nigsz=d->nvkat[i];
		ElemsForMats[ElemsIgMats[nigsz]+ElemsSizeForMats[nigsz]]=i;
		ElemsSizeForMats[nigsz]++;
	}

	logfile<<"All maps done"<<endl;

	ofstream ofp;

	ofp.open("urm",ios::binary);
	ofp<nx<ny<nz;
	for(i=0;i<nx;i++){ofp<d->Xcrd[i];}
	for(i=0;i<ny;i++){ofp<d->Ycrd[i];}
	for(i=0;i<nz;i++){ofp<d->Zcrd[i];}
	ofp.close();
	ofp.clear();

	ofp.open("urmm",ios::binary);
	for(i=0;i<kpar_reg;i++){ofp<d->nvkat[d->regular[i]];}
	ofp.close();
	ofp.clear();

	ofp.open("urms",ios::binary);
	for(lk=0;lk<nz;lk++)
	{
		for(lj=0;lj<ny;lj++)
		{
			for(li=0;li<nx;li++)
			{
				ofp<MatsSizeForNodes[nx*(ny*lk+lj)+li];
			}
		}
	}
	ofp.close();
	ofp.clear();

	logfile<<"Output regular mesh info done"<<endl;

	bool fnorm;
	int nn;
	int imat,mat,nelm,ielm,elm,elm_old,sz,ie[8],num_elem;
	float EpN[3][2], Ep[3][2], en[8][3];
	double invh, invh2, bf[8], x[8], y[8], z[8], rx[2], ry[2], rz[2], p[3], val, w;
	double ves_j[12][2], ves_j1[12][2], ves_j2[12][2];
	double val_j[3][2], val_j1[3][2], val_j2[3][2];
	double local_coords[3];
	nx=d->N_X;
	ny=d->N_Y;
	nz=d->N_Z;
	kx=nx-1;
	ky=ny-1;
	kz=nz-1;
	w=2.0*PI*d->nu;
	ofp.open("enormharm",ios::binary);
	for(lk=0;lk<nz;lk++)
	{
		p[2]=d->Zcrd[lk];
		for(lj=0;lj<ny;lj++)
		{
			p[1]=d->Ycrd[lj];
			for(li=0;li<nx;li++)
			{
				p[0]=d->Xcrd[li];				
				nn=nx*(ny*lk+lj)+li;
				nmat=MatsSizeForNodes[nn];
				nelm=ElemsSizeForNodes[nn];
				for(imat=0;imat<nmat;imat++)
				{
					Ep[0][0]=Ep[1][0]=Ep[2][0]=0.0;
					Ep[0][1]=Ep[1][1]=Ep[2][1]=0.0;
					mat=MatsForNodes[MatsIgNodes[nn]+imat];
					fnorm=true;
					elm_old=-1;
					num_elem=0;
					for(ielm=0;ielm<nelm;ielm++)
					{
						elm=ElemsForNodes[ElemsIgNodes[nn]+ielm];
						if(d->nvkat[elm]==mat)
						{
							if(!(d->nver[elm][13]>30))
							{
								if(elm!=elm_old)
								{
									for (i=0; i<12; i++)
									{
										e = tmap->ed[elm][i];
										ves_j[i][0]=u[2*e];
										ves_j[i][1]=u[2*e+1];
									}

									x[0]= d->xyz[d->nver[elm][0]][0];
									x[1]= d->xyz[d->nver[elm][7]][0];
									y[0]= d->xyz[d->nver[elm][0]][1];
									y[1]= d->xyz[d->nver[elm][7]][1];
									z[0]= d->xyz[d->nver[elm][0]][2];
									z[1]= d->xyz[d->nver[elm][7]][2];
									
									invh2=2.0/((x[1]-x[0])*(y[1]-y[0])*(z[1]-z[0]));
								}

								rx[0]=x[1]-p[0];
								rx[1]=p[0]-x[0];

								ry[0]=y[1]-p[1];
								ry[1]=p[1]-y[0];

								rz[0]=z[1]-p[2];
								rz[1]=p[2]-z[0];

								val_j[0][0]=(ves_j[0][0]*ry[0]*rz[0]+ves_j[1][0]*ry[1]*rz[0]+ves_j[2][0]*ry[0]*rz[1]+ves_j[3][0]*ry[1]*rz[1])*invh2;
								val_j[1][0]=(ves_j[4][0]*rx[0]*rz[0]+ves_j[5][0]*rx[0]*rz[1]+ves_j[6][0]*rx[1]*rz[0]+ves_j[7][0]*rx[1]*rz[1])*invh2;
								val_j[2][0]=(ves_j[8][0]*rx[0]*ry[0]+ves_j[9][0]*rx[1]*ry[0]+ves_j[10][0]*rx[0]*ry[1]+ves_j[11][0]*rx[1]*ry[1])*invh2;

								val_j[0][1]=(ves_j[0][1]*ry[0]*rz[0]+ves_j[1][1]*ry[1]*rz[0]+ves_j[2][1]*ry[0]*rz[1]+ves_j[3][1]*ry[1]*rz[1])*invh2;
								val_j[1][1]=(ves_j[4][1]*rx[0]*rz[0]+ves_j[5][1]*rx[0]*rz[1]+ves_j[6][1]*rx[1]*rz[0]+ves_j[7][1]*rx[1]*rz[1])*invh2;
								val_j[2][1]=(ves_j[8][1]*rx[0]*ry[0]+ves_j[9][1]*rx[1]*ry[0]+ves_j[10][1]*rx[0]*ry[1]+ves_j[11][1]*rx[1]*ry[1])*invh2;
							}
							else{
								if(elm!=elm_old)
								{
									for (i=0; i<12; i++)
									{
										e = tmap->ed[elm][i];
										ves_j2[i][0]=ves_j1[i][0]=ves_j[i][0]=u[2*e];
										ves_j2[i][1]=ves_j1[i][1]=ves_j[i][1]=u[2*e+1];
									}
									for (i=0; i<8; i++)
									{
										x[i] = d->xyz[d->nver[elm][i]][0];
										y[i] = d->xyz[d->nver[elm][i]][1];
										z[i] = d->xyz[d->nver[elm][i]][2];
									}
								}
								T_Brick L(x, y, z, d->nver[elm][13]);
					
								local_coords[0]=2.0*RegularLocalCoord[3*(ElemsIgNodes[nn]+ielm)]-1;
								local_coords[1]=2.0*RegularLocalCoord[3*(ElemsIgNodes[nn]+ielm)+1]-1;
								local_coords[2]=2.0*RegularLocalCoord[3*(ElemsIgNodes[nn]+ielm)+2]-1;

								L.Calc_value_inside_hex(ves_j[0], local_coords, val_j[0]);
								L.Calc_value_inside_hex(ves_j1[0], local_coords, val_j1[0]);
								L.Calc_value_inside_hex(ves_j2[0], local_coords, val_j2[0]);

								L.Calc_value_inside_hex(ves_j[1], local_coords, val_j[1]);
								L.Calc_value_inside_hex(ves_j1[1], local_coords, val_j1[1]);
								L.Calc_value_inside_hex(ves_j2[1], local_coords, val_j2[1]);
							}
							Ep[0][0] += w * val_j[0][1];
							Ep[1][0] += w * val_j[1][1];
							Ep[2][0] += w * val_j[2][1];
							Ep[0][1] += -w * val_j[0][0];
							Ep[1][1] += -w * val_j[1][0];
							Ep[2][1] += -w * val_j[2][0];
							num_elem++;

							// Добавить нормальное один раз
							if(fnorm)
							{
								for(lt=0;lt<2;lt++)
								{
									for(k=0;k<8;k++)
									{
										sz=(int)nodesLine.item[d->nver[elm][k]].size();
										for(j=0;j<sz;j++)
										{
											if(mat==nodesLine.item[d->nver[elm][k]][j].material)
											{
												ie[k]=j;
												en[k][0]=nodesLine.item[d->nver[elm][k]][j].val[0][lt];
												en[k][1]=nodesLine.item[d->nver[elm][k]][j].val[1][lt];
												en[k][2]=nodesLine.item[d->nver[elm][k]][j].val[2][lt];
												break;
											}
										}
										if(j==sz)break;
									}

									if(k==8)
									{
										if(!(d->nver[elm][13]>30))
										{
											if(elm!=elm_old)
											{
												x[0]= d->xyz[d->nver[elm][0]][0];
												x[1]= d->xyz[d->nver[elm][7]][0];
												y[0]= d->xyz[d->nver[elm][0]][1];
												y[1]= d->xyz[d->nver[elm][7]][1];
												z[0]= d->xyz[d->nver[elm][0]][2];
												z[1]= d->xyz[d->nver[elm][7]][2];

												invh=1.0/((x[1]-x[0])*(y[1]-y[0])*(z[1]-z[0]));
											}
												
											rx[0]=x[1]-p[0];
											rx[1]=p[0]-x[0];

											ry[0]=y[1]-p[1];
											ry[1]=p[1]-y[0];

											rz[0]=z[1]-p[2];
											rz[1]=p[2]-z[0];

											bf[0]=rx[0]*ry[0]*rz[0]*invh;
											bf[1]=rx[1]*ry[0]*rz[0]*invh;
											bf[2]=rx[0]*ry[1]*rz[0]*invh;
											bf[3]=rx[1]*ry[1]*rz[0]*invh;
											bf[4]=rx[0]*ry[0]*rz[1]*invh;
											bf[5]=rx[1]*ry[0]*rz[1]*invh;
											bf[6]=rx[0]*ry[1]*rz[1]*invh;
											bf[7]=rx[1]*ry[1]*rz[1]*invh;

											EpN[0][lt]=EpN[1][lt]=EpN[2][lt]=0.0;
											for(i=0;i<8;i++)
											{
												EpN[0][lt]+=bf[i]*en[i][0];
												EpN[1][lt]+=bf[i]*en[i][1];
												EpN[2][lt]+=bf[i]*en[i][2];
											}
										}
										else
										{
											local_coords[0]=RegularLocalCoord[3*(ElemsIgNodes[nn]+ielm)];
											local_coords[1]=RegularLocalCoord[3*(ElemsIgNodes[nn]+ielm)+1];
											local_coords[2]=RegularLocalCoord[3*(ElemsIgNodes[nn]+ielm)+2];
											for(i=0;i<8;i++)ves_j[i][lt]=en[i][0];
											GetNodeSolutionByLocalcCoords(local_coords,ves_j[lt],val);
											EpN[0][lt]=val;
											for(i=0;i<8;i++)ves_j[i][lt]=en[i][1];
											GetNodeSolutionByLocalcCoords(local_coords,ves_j[lt],val);
											EpN[1][lt]=val;
											for(i=0;i<8;i++)ves_j[i][lt]=en[i][2];
											GetNodeSolutionByLocalcCoords(local_coords,ves_j[lt],val);
											EpN[2][lt]=val;
										}
										fnorm=false;
									}
								}
							}
						}
						elm_old=elm;
					}
					Ep[0][0]/=num_elem;
					Ep[1][0]/=num_elem;
					Ep[2][0]/=num_elem;
					Ep[0][1]/=num_elem;
					Ep[1][1]/=num_elem;
					Ep[2][1]/=num_elem;
					ofp<mat;
					ofp<(Ep[0][0]+EpN[0][0])<(Ep[1][0]+EpN[1][0])<(Ep[2][0]+EpN[2][0]);
					ofp<(Ep[0][1]+EpN[0][1])<(Ep[1][1]+EpN[1][1])<(Ep[2][1]+EpN[2][1]);
				}
			}
		}
	}
	ofp.close();
	ofp.clear();
}
//------------------------------------------------------------------------
int main()
{
	int i;
	ifstream inf;
	int StartType,nproc;
	double tmpd;

	logfile.open("logharm3dOutput");

	GPDocSettings.Read("gpsettings");

	Vec_Prep_Data *d=NULL;
	if((d = new Vec_Prep_Data())==0) 
		throw logic_error("no memory for 3d mesh");

	inf.open("Harm3DParams");
	inf>>StartType;
	inf>>d->tasktype;
	inf.close();
	inf.clear();

	omp.InitNumThreads();
	nproc = (StartType==2)? 1 : omp.GetMaxThreads()-1;
	omp.SetNumThreads(nproc);

	// 0 - prosto raschet
	// 1 - s vikladkoi summarnogo polya po uzlam
	// 2 - dlya rascheta polya vliyaniya

	
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

	if(d->tasktype==0)
	{
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

	if(d->tasktype!=0)
	{
		if (d->LoadVectorE0ForLine(tmap->n)!=0)
			throw logic_error("Can't read normal field.");
	}

	double *u=NULL;

	if((u = new double[2*(tmap->n)])==0) Memory_allocation_error("u", "main");
	
	r.Read_Bin_File_Of_Double("v3.dat", u, tmap->n_c, 2);

	// выдача результата

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


	if(d->tasktype!=0 || StartType)
	{
		d->LoadReceivers("xyzVectorB", d->n_pointresB, d->pointresB);
		d->n_pointresE=0;
		d->n_pointres=d->n_pointresB;
		give_out = new Give_out_vec_mt(d, tmap, u);
		if(d->n_pointresB)
		{
			if(!GPDocSettings.vfem_direct_output)
			{
				d->xyzt=d->xyz;
				if(xyz0)
				{
					d->xyz=xyz0;
				}
				give_out->resultantB = new OutputResultant3d(give_out,vtWithoutDiscontinuity);
				give_out->resultantB->Prepare();
				if(xyz0)
				{
					d->xyz=d->xyzt;
					d->xyzt=NULL;
				}
			}	
			give_out->Give_out_on_hex(false);
			if(!GPDocSettings.vfem_direct_output)
			{
				give_out->resultantB->StopSolvers();
				delete give_out->resultantB;
				give_out->resultantB = NULL;
			}
		}
		give_out->Write_B_to_files_for_harm_loop(StartType,d->fdirect);
		delete give_out;

		d->LoadReceivers("xyzVectorE", d->n_pointresE, d->pointresE);
		d->n_pointresB=0;
		d->n_pointres=d->n_pointresE;
		give_out = new Give_out_vec_mt(d, tmap, u);
		if(d->n_pointresE)
		{
			if(!GPDocSettings.vfem_direct_output)
			{
				d->xyzt=d->xyz;
				if(xyz0)
				{
					d->xyz=xyz0;
				}
				give_out->resultantA = new OutputResultant3d(give_out,vtWithDiscontinuity);
				give_out->resultantA->Prepare();
				if(xyz0)
				{
					d->xyz=d->xyzt;
					d->xyzt=NULL;
				}
			}
			give_out->Give_out_on_hex(false);
			if(!GPDocSettings.vfem_direct_output)
			{
				give_out->resultantA->StopSolvers();
				delete give_out->resultantA;
				give_out->resultantA = NULL;
			}
		}
		give_out->Write_E_to_files_for_harm_loop(StartType,d->fdirect);
		delete give_out;
	}
	else
	{
		// МТЗ
		give_out = new Give_out_vec_mt(d, tmap, u);
		if(!GPDocSettings.vfem_direct_output)
		{
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
		}
		give_out->Give_out_on_hex(true);
		if(!GPDocSettings.vfem_direct_output)
		{
			give_out->resultantB->StopSolvers();
			delete give_out->resultantB;
			give_out->resultantB = NULL;
			give_out->resultantA->StopSolvers();
			delete give_out->resultantA;
			give_out->resultantA = NULL;
		}
		delete give_out;
	}

	if(StartType==1 && d->tasktype!=0)
	{
		OutEnormHarm(d,tmap,u);
	}

	if(xyz0){delete [] xyz0;xyz0=NULL;}

	if(u) {delete [] u; u=NULL;}
	if(d) {delete d; d=NULL;}
	if(tmap) {delete tmap; tmap=NULL;}

	logfile.close();

	return 0;
}
