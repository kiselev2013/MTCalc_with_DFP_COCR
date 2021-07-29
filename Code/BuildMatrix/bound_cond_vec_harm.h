#pragma once

class Bound_cond_vec_harm // учёт краевых условий
{
public:
	bool *isBlockBound;
	long *bound_nodes;
	long (*edges)[2];

	long n_edges;  
	long n_nodes; 
	long n_elem;
	long n_bound_nodes;

	long n_bound_edges;

	Bound_cond_vec_harm(long n_nodes, long n_edges, long n_bound_nodes,
		long *bound_nodes,  long (*edges)[2], long n_elem);
	~Bound_cond_vec_harm();

	void MakeListBoundsBlocks(); // алгоритм не совсем правильный
	void SetHomogenDirichletCond(long *ig, long *jg, long *idi, long *ijg,  
		double *di, double *gg, double *pr);

	//---------------AV formulation----------------------------------------
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
	Bound_cond_vec_harm(long n_nodes, long n_nodes_c, long n_edges, long n_edges_c, long n_bound_nodes,
		long *bound_nodes, long (*edges)[2], long *nded, long *nded_type, long *nvkat, 
		long (*nver)[14], long (*ed)[25], long *nodes_position_in_nded, long *edges_position_in_nded, long unk_c, long n_elem);
	void FormListM_bound();
	void SetHomogenDirichletCond_AV(long *ig, long *jg, long *idi, long *ijg,  
		double *di, double *gg, double *pr);
	//---------------AV formulation----------------------------------------
};