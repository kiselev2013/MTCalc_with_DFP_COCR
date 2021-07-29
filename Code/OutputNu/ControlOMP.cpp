#include "stdafx.h"
#include "ControlOMP.h"
//------------------------------------------------------------------------
ControlOMP::ControlOMP()
{
	isInit = false;
	NUMBER_OF_THREADS = 0;
	MAX_THREADS = 0;
	nMinDotProduct = 1000;
	nMinSparseMultMV = 10;
	nMinSparseMultMV2x2 = 10;
}
//------------------------------------------------------------------------
ControlOMP::~ControlOMP()
{
}
//------------------------------------------------------------------------
void ControlOMP::InitNumThreads()
{
	if (!isInit)
	{
		NUMBER_OF_THREADS = 1;
		MAX_THREADS = omp_get_max_threads();
		isInit = true;
	}
}
//------------------------------------------------------------------------
void ControlOMP::SetNumThreads(int num)
{
	if (isInit)
	{
		if (num <= MAX_THREADS)
			NUMBER_OF_THREADS = num;
	}
}
//------------------------------------------------------------------------
int ControlOMP::GetNumberOfThreads()
{
	return NUMBER_OF_THREADS;
}
//------------------------------------------------------------------------
int ControlOMP::GetMaxThreads()
{
	return MAX_THREADS;
}
//------------------------------------------------------------------------
int ControlOMP::GetNMinSparseMultMV()
{
	return nMinSparseMultMV;
}
//------------------------------------------------------------------------
int ControlOMP::GetNMinSparseMultMV2x2()
{
	return nMinSparseMultMV2x2;
}
//------------------------------------------------------------------------
int ControlOMP::GetNMinDotProduct()
{
	return nMinDotProduct;
}
//------------------------------------------------------------------------
