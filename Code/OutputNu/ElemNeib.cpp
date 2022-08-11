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
 *  Functions for reading and writing a structure of elements neighbors
 *
 *  Written by Prof. Marina G. Persova
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  mpersova@mail.ru
 *  Version 1.5 December 17, 2020
*/
#include "stdafx.h"
#include "ElemNeib.h"
//-----------------------------------------------------------
// Reading a structure of elements neighbors from file
//-----------------------------------------------------------
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
//-----------------------------------------------------------
// Writing a structure of elements neighbors to file
//-----------------------------------------------------------
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
