#pragma once
//------------------------------------------------------------------------
class FormatConverter
{
public:
	void FromRSFToCSR_Real_1_Sym(int nb, int *ig, int *sz_iptr, int *sz_jptr);
	void From2x2ToCSR_Complex_1_Sym(int nb, int *ig, int *idi, int *ijg,int *sz_iptr, int *sz_jptr);
	void From2x2ToCSR_Complex_1_NS(int nb, int *ig, int *idi, int *ijg,int *sz_iptr, int *sz_jptr);
		void From2x2ToCSRComplex_2_Sym(int nb, int *ig, int *jg, int *idi, int *ijg,double *di_block, double *ggl_block,
			MKL_INT64 *iptr, MKL_INT64 *jptr, double *aelem);
	void FromRSFToCSR_Real_2_Sym(int nb, int *ig, int *jg, double *di, double *gg,MKL_INT64 *iptr, MKL_INT64 *jptr, double *aelem);
	void From2x2ToCSRComplex_2_NS(int nb, int *ig, int *jg, int *idi, int *ijg,double *di_block, double *ggl_block, double *ggu_block,
		MKL_INT64 *iptr, MKL_INT64 *jptr, double *aelem);
	FormatConverter();
	~FormatConverter();
};
