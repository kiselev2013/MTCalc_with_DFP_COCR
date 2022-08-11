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
 *  The file bound_cond_vec_harm.h contains the headers for the functions of the class "Bound_cond_vec_harm"        
 *  for imposing the boundary conditions in the 3D vector FEM time-harmonic global matrix and right hand side vector
 *  in the case of both single (A) and coupled potentials (AV)                                                      
 *                                                                                                                  
 *  Written by Ph.D. Petr A. Domnikov                                                                               
 *  Novosibirsk State Technical University,                                                                         
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                               
 *  p_domnikov@mail.ru                                                                                              
 *  Version 1.2 January 11, 2021                                                                                    
*/                                                                                                                  


#pragma once

//-----------------------------------------------------------------------------
// The class for imposing the boundary conditions in the 3D vector FEM time-harmonic global matrix and right hand side vector
//-----------------------------------------------------------------------------
class Bound_cond_vec_harm 
{
public:
	bool *isBlockBound; // the list of flags if the boundary condition is imposed to the edge 
	long *bound_nodes;  // the list of boundary nodes
	long (*edges)[2];   // the list of edges by two verteces

	long n_edges;        // the number of edges in the mesh
	long n_nodes;        // the number of nodes in the mesh
	long n_elem;         // the number of finite elements in the mesh
	long n_bound_nodes;  // the number of boundary nodes in the mesh
	long n_bound_edges;  // the number of boundary edges in the mesh

	// Constructor for the class Bound_cond_vec_harm (single potential A)
	Bound_cond_vec_harm(long n_nodes, long n_edges, long n_bound_nodes,
		long *bound_nodes,  long (*edges)[2], long n_elem);

	// Destructor for the class Bound_cond_vec_harm
	~Bound_cond_vec_harm();

	// Setting the flags if the boundary condition is imposed to nodes and edges
	void MakeListBoundsBlocks(); 

	// Imposing the homogenious Dirichlet boundary conditions for the single potential A
	void SetHomogenDirichletCond(long *ig, long *jg, long *idi, long *ijg,  
		double *di, double *gg, double *pr);

	// variables and functions for the coupled potentials A-V

	bool *m_bound;
	long *nded;
	long *nded_type;
	long *nvkat;
	long (*nver)[14];
	long (*ed)[25];
	long *nodes_position_in_nded;
	long *edges_position_in_nded;
	long n_nodes_c;
	long n_edges_c;
	long unk_c;

	// Constructor for the class Bound_cond_vec_harm (coupled potentials A-V)
	Bound_cond_vec_harm(long n_nodes, long n_nodes_c, long n_edges, long n_edges_c, long n_bound_nodes,
		long *bound_nodes, long (*edges)[2], long *nded, long *nded_type, long *nvkat, 
		long (*nver)[14], long (*ed)[25], long *nodes_position_in_nded, long *edges_position_in_nded, long unk_c, long n_elem);
	void FormListM_bound();

	// Imposing the homogenious Dirichlet boundary conditions for the coupled potentials AV
	void SetHomogenDirichletCond_AV(long *ig, long *jg, long *idi, long *ijg,  
		double *di, double *gg, double *pr);
};
