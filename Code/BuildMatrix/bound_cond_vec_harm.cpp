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
 * The file bound_cond_vec_harm.cpp contains functions of the class "Bound_cond_vec_harm"                            
 *  for imposing the boundary conditions in the 3D vector FEM time-harmonic global matrix and right hand side vector  
 *  in the case of both single (A) and coupled potentials (AV)                                                                                         
 *                                         
 *  Written by Ph.D. Petr A. Domnikov                                                     
 *  Novosibirsk State Technical University,                                               
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                     
 *  p_domnikov@mail.ru                                                                    
 *  Version 1.2 January 11, 2021                                                          
*/                                                                                        


#include "stdafx.h"
#include "bound_cond_vec_harm.h"
//-----------------------------------------------------------------------------
// Constructor for the class Bound_cond_vec_harm (single potentials A-V)
//-----------------------------------------------------------------------------
Bound_cond_vec_harm::Bound_cond_vec_harm(long n_nodes, long n_edges, long n_bound_nodes,
						   long *bound_nodes,  long (*edges)[2], long n_elem)
{
	this->n_nodes = n_nodes;
	this->n_edges = n_edges;
	this->n_bound_nodes = n_bound_nodes;
	this->bound_nodes = bound_nodes;
	this->edges = edges;
	this->n_elem = n_elem;

	isBlockBound = m_bound = NULL;

	MakeListBoundsBlocks();
}
//-----------------------------------------------------------------------------
// Constructor for the class Bound_cond_vec_harm (coupled potentials A-V)
//-----------------------------------------------------------------------------
Bound_cond_vec_harm::Bound_cond_vec_harm(long n_nodes, long n_nodes_c, long n_edges, long n_edges_c, long n_bound_nodes,
										 long *bound_nodes, long (*edges)[2], long *nded, long *nded_type, long *nvkat, 
										 long (*nver)[14], long (*ed)[25], long *nodes_position_in_nded, long *edges_position_in_nded, long unk_c, long n_elem)
{
	this->n_nodes = n_nodes;
	this->n_nodes_c = n_nodes_c;
	this->n_edges = n_edges;
	this->n_edges_c = n_edges_c;
	this->n_bound_nodes = n_bound_nodes;
	this->bound_nodes = bound_nodes;
	this->edges = edges;
	this->nded = nded;
	this->nded_type = nded_type;
	this->nvkat = nvkat;
	this->nver = nver;
	this->ed = ed;
	this->nodes_position_in_nded = nodes_position_in_nded;
	this->edges_position_in_nded = edges_position_in_nded;
	this->n_elem = n_elem;

	this->unk_c = unk_c;
	isBlockBound = m_bound = NULL;

	MakeListBoundsBlocks();
}
//-----------------------------------------------------------------------------
// Desstructor for the class Bound_cond_vec_harm
//-----------------------------------------------------------------------------
Bound_cond_vec_harm::~Bound_cond_vec_harm()
{
	if(isBlockBound) {delete [] isBlockBound; isBlockBound=NULL;}
	if(m_bound) {delete [] m_bound; m_bound=NULL;}
}
//------------------------------------------------------------------------
// Setting the flags if the boundary condition is imposed to nodes and edges
//------------------------------------------------------------------------
void Bound_cond_vec_harm::MakeListBoundsBlocks()
{
	long i;
	bool *is_node_bound=NULL;

	// fill the is_node_bound[]
	if((is_node_bound = new bool[n_nodes])==0) 
		Memory_allocation_error("is_node_bound","Bound_cond_vec_harm::Make_list_of_bound_edges_harm");

	for(i=0; i<n_nodes; i++)
		is_node_bound[i] = false;

	for(i=0; i<n_bound_nodes; i++)
		is_node_bound[bound_nodes[i]] = true;

	// fill the isBlockBound[]
	if((isBlockBound = new bool[n_edges])==0)
		Memory_allocation_error("isBlockBound","Bound_cond_vec_harm::Make_list_of_bound_edges_harm");;

	for(i=0; i<n_edges; i++)
		isBlockBound[i] = false;

	for(i=0; i<n_edges; i++)
		if(is_node_bound[edges[i][0]]==1 && is_node_bound[edges[i][1]]==1)		
			isBlockBound[i] = true;

	if(is_node_bound) {delete [] is_node_bound; is_node_bound=NULL;}
}
//------------------------------------------------------------------------
void Bound_cond_vec_harm::FormListM_bound()
{
	m_bound = new bool[unk_c];

	// calculating the number of boundary edges
	n_bound_edges = 0;
	for(long i = 0; i < n_edges; i++)
		if(isBlockBound[i]) n_bound_edges++;

	for(long ik = 0; ik < unk_c; ik++) m_bound[ik] = true;

	for(long ielem = 0; ielem < n_elem; ielem++)
	{
		if(true)
		{
			for(int ii = 0; ii < 8; ii++)
			{
				if(nver[ielem][ii] < n_nodes_c)
					m_bound[nodes_position_in_nded[nver[ielem][ii]]] = false;
			}
			for(int ii = 0; ii < 12; ii++)
			{
				if(ed[ielem][ii] < n_edges_c)
					m_bound[edges_position_in_nded[ed[ielem][ii]]] = false;
			}
		}
		else
		{
			for(int ii = 0; ii < 12; ii++)
			{
				if(ed[ielem][ii] < n_edges_c)
					m_bound[edges_position_in_nded[ed[ielem][ii]]] = false;
			}
		}
	}

	for(long iedg = 0; iedg < n_edges; iedg++)
	{
		if(isBlockBound[iedg])
		{
			if(iedg < n_edges_c)
			{
				for(long i = 0; i < unk_c; i++)
				{
					if(nded[i] == iedg)
					{
						if(nded_type[i] == 1)
						{
							m_bound[i] = true;
							break;
						}
					}
				}
			}
		}
	}

	for(long inode = 0; inode < n_bound_nodes; inode++)
	{
		if(bound_nodes[inode] < n_nodes_c)
		{
			for(long i = 0; i < unk_c; i++)
			{
				if(nded[i] == bound_nodes[inode])
				{
					if(nded_type[i] == 0)
					{
						m_bound[i] = true;
						break;
					}
				}
			}
		}
	}
}
//------------------------------------------------------------------------
// Imposing the homogenious Dirichlet boundary conditions for the coupled potentials A-V
//------------------------------------------------------------------------
void Bound_cond_vec_harm::SetHomogenDirichletCond_AV(long *ig, long *jg, long *idi, long *ijg,  
													 double *di, double *gg, double *pr)
{
	long i, j, k, adr;
	long size;
	bool flag;

	cout << "Set_dirichlet_cond_harm AV...\n";

	long nded_size = n_edges + n_nodes;

	for(i = 0; i < nded_size; i++)
	{
		flag = false;
		if(nded_type[i] == 1 && nded[i] < n_edges_c) flag = true;
		else if(nded_type[i] == 0 && nded[i] < n_nodes_c) flag = true;

		if(flag)
		{
			if(m_bound[i])
			{
				di[idi[i]] = 1.0;
				size = idi[i+1] - idi[i];
				if(size==2)
					di[idi[i]+1] = 0.0;

				//gg
				for (j=ig[i]; j<=ig[i+1]-1; j++)
				{
					k = jg[j];
					adr = ijg[j];
					gg[adr] = 0.0;
					size = ijg[j+1] - adr;
					if(size==2)
						gg[adr+1] = 0.0;
				}

				//pr
				pr[i*2] = pr[i*2+1] = 0.0;
			}
		}
	}

	for(i = 0; i < nded_size; i++)
	{
		flag = false;
		if(nded_type[i] == 1 && nded[i] < n_edges_c) flag = true;
		else if(nded_type[i] == 0 && nded[i] < n_nodes_c) flag = true;

		if(flag)
		{
			for(j = ig[i]; j < ig[i + 1]; j++)
			{
				k = jg[j];
				if(m_bound[k] == true)
				{
					if(m_bound[i] == false)
					{
						adr = ijg[j];
						gg[adr] = 0.0;
						size = ijg[j+1] - adr;
						if(size==2)
							gg[adr+1] = 0.0;
					}
				}
			}
		}
	}
}
//------------------------------------------------------------------------
// Imposing the homogenious Dirichlet boundary conditions for the single potential A
//------------------------------------------------------------------------
void Bound_cond_vec_harm::SetHomogenDirichletCond(long *ig, long *jg, long *idi, long *ijg,  
	double *di, double *gg, double *pr)
{
	long i, j, k, adr;
	long size;

	cout << "Set_dirichlet_cond_harm...\n";

	for (i=0; i<n_edges; i++)
	{
		if (isBlockBound[i])
		{
			// diagonal elements processing
			di[idi[i]] = 1.0;
			size = idi[i+1] - idi[i];
			if(size==2)
				di[idi[i]+1] = 0.0;

			// non-diagonal elements processing
			for (j=ig[i]; j<=ig[i+1]-1; j++)
			{
				k = jg[j];
				adr = ijg[j];
				gg[adr] = 0.0;
				size = ijg[j+1] - adr;
				if(size==2)
					gg[adr+1] = 0.0;
			}

			// right hand side
			pr[i*2] = pr[i*2+1] = 0.0;
		}// end if

		// symmetrization
		for (j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			if (isBlockBound[k])
			{
				adr = ijg[j];
				gg[adr] = 0.0;
				size = ijg[j+1] - adr;
				if(size==2)
					gg[adr+1] = 0.0;
			}
		}// end j
	}// end i
}
//------------------------------------------------------------------------