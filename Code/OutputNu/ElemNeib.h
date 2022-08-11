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
 *  Header file for ElemNeib.cpp
 *
 *  Written by Prof. Marina G. Persova
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  mpersova@mail.ru
 *  Version 1.5 December 17, 2020
*/
#pragma once

struct ElemNeib
{
	vector<int> neib[6]; // lists of neighbors for element surfaces
	ElemNeib()
	{
		neib[0].clear();neib[1].clear();neib[2].clear();
		neib[3].clear();neib[4].clear();neib[5].clear();
	}
};	

int ReadElemNeib(vector<ElemNeib> &ElemNeibVec,int kpar);
int WriteElemNeib(vector<ElemNeib> &ElemNeibVec,int kpar);
