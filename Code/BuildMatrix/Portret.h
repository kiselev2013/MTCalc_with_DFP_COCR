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
 *  This file contains some basic routines for construction of a matrix portrait of a finite element SLAE 
 * (nodal FEM for solving 2D and 3D problems) 
 *                                                                                               
 *  Written by Ph.D. Petr A. Domnikov                                                          
 *  Novosibirsk State Technical University,                                                    
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                          
 *  p_domnikov@mail.ru                                                                         
 *  Version 1.3 April 5, 2021                                                                  
*/                                                                                             
#pragma once
//----------------------------
struct list // list of integers
{
	long number;
	list *next;
};
//----------------------------
class Portret
{
public:

	long *ig;
	long *jg;
	list *s; // structure
	bool is_mem_ig_jg_allocated;
	long n; // dimension of the global matrix
	long n_of_elements;
	long size_jg;
	long (*nvtr)[8];
	long (*nvtr4)[4];
	long (*nver)[14];

	long n_c; // number of non-terminal nodes
	long *ig_t; // T-transform matrix
	long *jg_t;
	
	Portret(long (*nver)[14], long n_of_elements, long n_of_nodes, long n_c, long *ig_t, long *jg_t);
	Portret(long (*nvtr)[4],  long n_of_elements, long n_of_nodes, long n_c, long *ig_t, long *jg_t); // for 2d
	Portret(long (*nvtr)[8],  long n_of_elements, long n_of_nodes, long n_c, long *ig_t, long *jg_t);
	Portret(long (*nvtr4)[4], long n_of_elements, long n_of_nodes);
	
	Portret();
	~Portret();

	void Gen_T_Portrait(); // for 3D
	void Gen_T_Portrait2(); // for 3D
	void Gen_Portret_2d_Linear(); // for spline
	void Gen_Portret_2d_rect_T(); // for 2D

	void Add_In_Ordered_List(list *s, long x); // inserting an element into an ordered list
	void Clear_List(list *s);
	long Size_Of_List(list *s);
};
