#pragma once

class T_Portrait
{
public:
	long (*ed)[25];
	long n_elem,n,size_jg;
	long *ig,*jg,*idi,*ijg; 
	T_Portrait(long *ed,long n,long n_elem);
	~T_Portrait();
	void Gen_Portrait();
	void Gen_idi_ijg(long *nvkat, long (*nver)[14]);
	void Set_type_of_block(long *target_array, long adr, long type);
};

const int FILTER_MASS_MATRIX_VEC[12][12] = { 
// 1  2  3  4  5  6  7  8  9  10 11 12
   2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
   2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
   2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
   2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 
   1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 
   1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 
   1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,  
   1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,  
   1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,  
   1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2
};
