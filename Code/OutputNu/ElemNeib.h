#pragma once


struct ElemNeib{
	vector<int> neib[6];
	ElemNeib(){
		neib[0].clear();neib[1].clear();neib[2].clear();
		neib[3].clear();neib[4].clear();neib[5].clear();
	}
};	


int ReadElemNeib(vector<ElemNeib> &ElemNeibVec,int kpar);
int WriteElemNeib(vector<ElemNeib> &ElemNeibVec,int kpar);

