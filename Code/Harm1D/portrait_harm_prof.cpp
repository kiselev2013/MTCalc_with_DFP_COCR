#include "stdafx.h"
#include "portrait_harm_prof.h"
//-----------------------------------------------------------
Portrait_profil::Portrait_profil(long n_elem)
{
	this->n_elem = n_elem;
	this->mem = false;
}
//-----------------------------------------------------------
Portrait_profil::~Portrait_profil()
{
	if(this->mem == true)
	{
		delete [] this->ig;
	}
}
//-----------------------------------------------------------
void Portrait_profil::Gen_ig_for_1d_harm_problem()
{
	long i, k;

	this->n = (this->n_elem + 1)*2;

	this->ig = new long[n + 1];
	if(ig == 0)
	{
		char var[] = {"ig"};
		char func[] = {"Portrait_profil::Gen_ig_for_1d_harm_problem"};
		printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
		exit(1);
	}
	this->mem = true;

	ig[0] = 0;
	ig[1] = 0;
	ig[2] = 1;
	k = 3;

	for(i=0; i<this->n_elem; i++)
	{
		ig[k] = ig[k-1] + 2;
		ig[k+1] = ig[k] + 3;
		k += 2; 		
	}
}
//-----------------------------------------------------------