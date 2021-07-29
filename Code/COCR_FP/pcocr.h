#pragma once
//------------------------------------------------------------------------
class PCOCR
{
public:
	int PCOCR_2x2_Folded(int n, int *ig, int *jg, int *idi, int *ijg, double *di, double *gg, double *pr,
		double *x, double eps, int maxiter, double *y_omp,
		int kpar, int n_edges_c, int n_nodes_c,
		double (*xyz)[3], int *nvkat, int (*nver)[14], double *sigma3d, int (*edges)[2],
		int *ig_t, int *jg_t, double *gg_t, int *is_node_bound);
	PCOCR();
	~PCOCR();
};
