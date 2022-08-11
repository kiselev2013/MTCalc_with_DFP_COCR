/**                                                                                                    
 * GENERAL REMARKS                                                                                     
 *                                                                                                     
 *  This code is freely available under the following conditions:                                      
 *                                                                                                     
 *  1) The code is to be used only for non-commercial purposes.                                        
 *  2) No changes and modifications to the code without prior permission of the developer.             
 *  3) No forwarding the code to a third party without prior permission of the developer.              
 *                                                                                                     
 *  			MTCalc_with_DFP_COCR                                                           
 *  This file contains subroutines for managing the number of OMP threads                              
 *                                                                                                     
 *  Written by Ph.D. Petr A. Domnikov                                                                  
 *  Novosibirsk State Technical University,                                                            
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                  
 *  p_domnikov@mail.ru                                                                                 
 *  Version 1.3 April 12, 2021                                                                         
*/                                                                                                     

#include "stdafx.h"
#include "ControlOMP.h"
//------------------------------------------------------------------------   
// Constructor                                                               
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
// Destructor                                                                  
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
