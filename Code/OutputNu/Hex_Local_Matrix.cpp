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
 *  Computation of local matrices on hexahedrons                                          
 *                                                                                        
 *  Written by Ph.D. Petr A. Domnikov                                                     
 *  Novosibirsk State Technical University,                                               
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                     
 *  p_domnikov@mail.ru                                                                    
 *  Version 1.2 March 10, 2021                                                            
*/                                                                                                                                                                                 

#include "stdafx.h"
#include "Hex_Local_Matrix.h"
#include "gauss_3.h"

extern void Mult_Plot(double *a, double *x, double *y, int n);
extern double Scal(const int &n, const double *v1, const double *v2);

extern ofstream logfile;
//---------------------------------------------------
Hex_Local_Matrix::~Hex_Local_Matrix()
{
}
//---------------------------------------------------
Hex_Local_Matrix::Hex_Local_Matrix(int num, long (*nver)[14], double (*xyz)[3])
{
	int i;

	for(i=0; i<8; i++)
	{
		x[i] = xyz[nver[num][i]][0];
		y[i] = xyz[nver[num][i]][1];
		z[i] = xyz[nver[num][i]][2];
	}

	number_of_element = num;
	type_of_hex = nver[num][13];
	JforParCalc = false;
}
//---------------------------------------------------
void Hex_Local_Matrix::Calc_J(int n_of_point)
{
	int i, j;

	if (type_of_hex <= 30) // parallelepiped
	{
		if(JforParCalc)
			return;

		hx = x[7] - x[0];
		hy = y[7] - y[0];
		hz = z[7] - z[0];

		J[0][0] = hx/2;
		J[1][1] = hy/2;
		J[2][2] = hz/2;
		J[0][1] = J[1][0] = J[0][2] = J[2][0] = J[1][2] = J[2][1] = 0.0;

		det_J = hx*hy*hz/8.0;
		det_J_abs = fabs(det_J);

		J_1_T[0][0] = 2.0/hx/det_J;
		J_1_T[1][1] = 2.0/hy/det_J;
		J_1_T[2][2] = 2.0/hz/det_J;
		J_1_T[0][1]=J_1_T[1][0]=J_1_T[0][2]=J_1_T[2][0]=J_1_T[1][2]=J_1_T[2][1]=0.0;

		JforParCalc = true;
	} 
	else // hexahedron
	{
	for(i=0; i<3; i++)
	for(j=0; j<3; j++)
		J[i][j] = 0.0;

	// elements of the Jacobian matrix
	for(i=0; i<8; i++)
	{
		J[0][0] += x[i]*gauss_3_d_phi[n_of_point][i][0]; //  d_xi[i];
		J[0][1] += x[i]*gauss_3_d_phi[n_of_point][i][1]; //  d_eta[i];
		J[0][2] += x[i]*gauss_3_d_phi[n_of_point][i][2]; //  d_zeta[i];
	
		J[1][0] += y[i]*gauss_3_d_phi[n_of_point][i][0]; //  d_xi[i];
		J[1][1] += y[i]*gauss_3_d_phi[n_of_point][i][1]; //  d_eta[i];
		J[1][2] += y[i]*gauss_3_d_phi[n_of_point][i][2]; //  d_zeta[i];
	
		J[2][0] += z[i]*gauss_3_d_phi[n_of_point][i][0]; //  d_xi[i];
		J[2][1] += z[i]*gauss_3_d_phi[n_of_point][i][1]; //  d_eta[i];
		J[2][2] += z[i]*gauss_3_d_phi[n_of_point][i][2]; //  d_zeta[i];
	}

	// calculate Jacobian (determinant)
	this->det_J = J[0][0]*J[1][1]*J[2][2] - J[0][0]*J[1][2]*J[2][1] + J[1][0]*J[2][1]*J[0][2]
	    - J[1][0]*J[0][1]*J[2][2] + J[2][0]*J[0][1]*J[1][2] - J[2][0]*J[1][1]*J[0][2];

	// Jacobian modulus
	this->det_J_abs = fabs(det_J);

	// matrix inverse to Jacobi matrix (and transposed)
	J_1_T[0][0] = (J[1][1]*J[2][2]-J[2][1]*J[1][2])/det_J;
	J_1_T[1][0] =  (J[2][1]*J[0][2]-J[0][1]*J[2][2])/det_J;
    J_1_T[2][0] = (J[0][1]*J[1][2]-J[1][1]*J[0][2])/det_J;
    J_1_T[0][1] = (-J[1][0]*J[2][2]+J[2][0]*J[1][2])/det_J;
    J_1_T[1][1] = (J[0][0]*J[2][2]-J[2][0]*J[0][2])/det_J;
    J_1_T[2][1] = (-J[0][0]*J[1][2]+J[1][0]*J[0][2])/det_J;
    J_1_T[0][2] = (J[1][0]*J[2][1]-J[2][0]*J[1][1])/det_J;
    J_1_T[1][2] = (-J[0][0]*J[2][1]+J[2][0]*J[0][1])/det_J;
    J_1_T[2][2] = (J[0][0]*J[1][1]-J[1][0]*J[0][1])/det_J;

	}
}
//---------------------------------------------------
void Hex_Local_Matrix::Calc_local_matrix_b_for_parallelepiped()
{
	hx = x[7] - x[0];
	hy = y[7] - y[0];
	hz = z[7] - z[0];

	double t2, t3, t4, t5, t6;

	t2 = hx*hy*hz;
	t3 = t2/27.0;
	t4 = t2/54.0;
	t5 = t2/108.0;
	t6 = t2/216.0;

	b[0][0] = t3;
	b[0][1] = t4;
	b[0][2] = t4;
	b[0][3] = t5;
	b[0][4] = t4;
	b[0][5] = t5;
	b[0][6] = t5;
	b[0][7] = t6;
	b[1][0] = t4;
	b[1][1] = t3;
	b[1][2] = t5;
	b[1][3] = t4;
	b[1][4] = t5;
	b[1][5] = t4;
	b[1][6] = t6;
	b[1][7] = t5;
	b[2][0] = t4;
	b[2][1] = t5;
	b[2][2] = t3;
	b[2][3] = t4;
	b[2][4] = t5;
	b[2][5] = t6;
	b[2][6] = t4;
	b[2][7] = t5;
	b[3][0] = t5;
	b[3][1] = t4;
	b[3][2] = t4;
	b[3][3] = t3;
	b[3][4] = t6;
	b[3][5] = t5;
	b[3][6] = t5;
	b[3][7] = t4;
	b[4][0] = t4;
	b[4][1] = t5;
	b[4][2] = t5;
	b[4][3] = t6;
	b[4][4] = t3;
	b[4][5] = t4;
	b[4][6] = t4;
	b[4][7] = t5;
	b[5][0] = t5;
	b[5][1] = t4;
	b[5][2] = t6;
	b[5][3] = t5;
	b[5][4] = t4;
	b[5][5] = t3;
	b[5][6] = t5;
	b[5][7] = t4;
	b[6][0] = t5;
	b[6][1] = t6;
	b[6][2] = t4;
	b[6][3] = t5;
	b[6][4] = t4;
	b[6][5] = t5;
	b[6][6] = t3;
	b[6][7] = t4;
	b[7][0] = t6;
	b[7][1] = t5;
	b[7][2] = t5;
	b[7][3] = t4;
	b[7][4] = t5;
	b[7][5] = t4;
	b[7][6] = t4;
	b[7][7] = t3;
}
//------------------------------------------------------------------------
void Hex_Local_Matrix::CalcMassMatrix()
{
	if (type_of_hex <= 30)
	{
		Calc_local_matrix_b_for_parallelepiped();
	} 
	else
	{
		int i, j, i1, j1;
		double gauss_3_mult;

		for(i=0; i<8; i++) // first set 0
		{
			for(j=0; j<8; j++)
				b[i][j] = 0.0;
		}

		for(i=0; i<27; i++) // by the number of integration points
		{
			Calc_J(i); // calculate Jacobian
			gauss_3_mult = gauss_3_A_all[i]*det_J_abs; // A_i*A_j*A_k*|J|

			for(i1=0; i1<8; i1++)
			{
				for(j1=0; j1<8; j1++)
					b[i1][j1] += gauss_3_phi[i][i1]*gauss_3_phi[i][j1]*gauss_3_mult;
			}//i1
		}// i
	}
}
//------------------------------------------------------------------------
