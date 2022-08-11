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
 *  This file contains some basic routines for output                                         
 *                                                                                        
 *  Written by Ph.D. Petr A. Domnikov                                                     
 *  Novosibirsk State Technical University,                                               
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                     
 *  p_domnikov@mail.ru                                                                    
 *  Version 1.2 December 10, 2020                                                            
*/                                                                                        


#pragma once

struct NodeForLine
{
	long material;
	double val[3][2]; // for a node, all 3 components are used, for an edge, only one
};
//------------------------------------------------------------------------
class ListOfNodesForLine
{
public:
	long kuzlov;
	long size; // number of nodes with repetitions 
	vector< vector<NodeForLine> > item;

	ListOfNodesForLine();
	~ListOfNodesForLine();
	void Init(long kuzlov);
	void Add(long node, long mat);

	void CalcSize();
};
