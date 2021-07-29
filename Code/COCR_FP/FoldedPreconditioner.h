#pragma once
//------------------------------------------------------------------------
class FoldedPreconditioner
{
public:
	int n; // размерность исходной системы
	int m; // размерность пространства

	double zero;

	int *iptr_g;
	int *jptr_g;
	double *gg_g;

	int *iptr_gt;
	int *jptr_gt;
	double *gg_gt;

	double *diag_a;
	double *diag_gagt;
	double *help;
	double *help2;
	double *help3;

	FoldedPreconditioner();
	~FoldedPreconditioner();

	void BuildGMatrix(int kpar, int n_edges_c, int n_nodes_c,
		int *nvkat, int (*nver)[14], double *sigma3d, int (*edges)[2],
		int *ig_t, int *jg_t, double *gg_t, int *is_node_bound);

	void Prepare(int kpar, int n_edges_c, int n_nodes_c,
		int *nvkat, int (*nver)[14], double *sigma3d, int (*edges)[2],
		int *ig_t, int *jg_t, double *gg_t, int *is_node_bound,
		int *ig, int *jg, double *di, double *gg);

	void BuildDiag(int *ig, int *jg, double *di, double *gg);
	void BuildDiag(double *di, double c);
	void ApplyPreconditioner(double *x, double *y);

	void BuildDiagComplex(int *ig, int *jg, int *idi, int *ijg, double *di, double *gg);
	void ApplyPreconditionerComplex(double *x, double *y);

	void FromRSFtoCSR_1(int n, int *ig, int *sz_iptr, int *sz_jptr);
	void FromRSFtoCSR_2(int n, int *ig, int *jg, double *di, double *ggl, double *ggu,
		int *iptr, int *jptr, double *aelem);

	void From2x2ToCSR2x2_1(int n, int *ig, int *idi, int *ijg,
		int *sz_iptr, int *sz_jptr, int *sz_ijptr, int *sz_aelem);
	void From2x2ToCSR2x2_2(int n, int *ig, int *jg, int *idi, int *ijg,
		double *di_block, double *gg_block, 
		int *iptr, int *jptr, int *ijptr, double *aelem);
};