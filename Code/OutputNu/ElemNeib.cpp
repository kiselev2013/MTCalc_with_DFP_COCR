#include "stdafx.h"
#include "ElemNeib.h"

int ReadElemNeib(vector<ElemNeib> &ElemNeibVec,int kpar)
{
	int i,j,k,st;
	ifstream inf;
	inf.open("elem_neib",ios::binary);
	if(!inf)return 1;
	ElemNeibVec.clear();
	ElemNeibVec.resize(kpar);
	for(i=0;i<kpar;i++){
		for(j=0;j<6;j++){
			inf>st;
			ElemNeibVec[i].neib[j].resize(st);
			for(k=0;k<st;k++)
			{
				inf>ElemNeibVec[i].neib[j][k];
				ElemNeibVec[i].neib[j][k]--;
			}
		}
	}
	inf.close();
	inf.clear();

	return 0;
}

int WriteElemNeib(vector<ElemNeib> &ElemNeibVec,int kpar)
{
	int i,j,k,m,st;
	ofstream ofp;

	ofp.open("elem_neib",ios::binary);
	if(!ofp)return 1;
	for(i=0;i<kpar;i++){
		for(j=0;j<6;j++){
			st=(int)ElemNeibVec[i].neib[j].size();
			ofp<st;
			for(k=0;k<st;k++)
			{
				m=ElemNeibVec[i].neib[j][k]+1;
				ofp<m;
			}
		}
	}
	ofp.close();
	ofp.clear();

	return 0;
}
