#include "stdafx.h"
#include "T_Mapping.h"
//----------------------------------------------------------------------------------------- 
T_Mapping_Vec::T_Mapping_Vec(int (*nver)[14], double (*xyz)[3], int kuzlov, int kpar)
{
	this->nver = nver;
	this->xyz = xyz;
	this->kuzlov = kuzlov;
	this->kpar = kpar;

	ed = NULL;
	edges = NULL;
	gg_t = NULL;
	ig_t = NULL;
	jg_t = NULL;
	ig_s = NULL;
	jg_s = NULL;
	s_val = NULL;
}
//-----------------------------------------------------------------------------------------
T_Mapping_Vec::~T_Mapping_Vec()
{
	if(ed) {delete [] ed; ed = NULL;}
	if(gg_t) {delete [] gg_t; gg_t = NULL;}
	if(ig_t) {delete [] ig_t; ig_t = NULL;}
	if(jg_t) {delete [] jg_t; jg_t = NULL;}
	if(ig_s) {delete [] ig_s; ig_s = NULL;}
	if(jg_s) {delete [] jg_s; jg_s = NULL;}
	if(s_val) {delete [] s_val; s_val = NULL;}
	if(edges) {delete [] edges; edges = NULL;}
}
//-----------------------------------------------------------------------
