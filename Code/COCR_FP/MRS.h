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
 *  Minimal residual smoothing (MRS)                                                        
 *                                                                                          
 *  Written by Ph.D. Petr A. Domnikov                                                       
 *  Novosibirsk State Technical University,                                                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                       
 *  p_domnikov@mail.ru                                                                      
 *  Version 1.2 April 7, 2021                                                               
*/                                                                                          



#pragma once
class MRS
{
private:
	int n;
	int nb;
	double *y;
	double *s;
	double *u;
	double w;
	double w_re;
	double w_im;
	double residual;
	bool isInit;
	double eps;

public:
	// Real-valued MRS
	void Mrs(double *x, double *r);

	// Complex-valued MRS
	void MrsCmplx(double *x, double *r);

	double* GetResidual();
	void GetResidual(double *r);
	double GetResidualNorm();
	void GetSolution(double *a);
	double GetW();
	double GetWRe();
	double GetWIm();

	void CalcResidual();

	MRS(int n);
	~MRS();
};


