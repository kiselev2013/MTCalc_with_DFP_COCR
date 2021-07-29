#include "stdafx.h"
#include "T_Brick.h"
#include "gauss_3.h"
#include "gauss_3_vec.h"
#include "gauss3vec2d.h"

extern ofstream logfile;

//-------------------------------------------------------------------------------
//--- Этот конструктор используется в том случае, если требуется только 
//--- выдать значение внутри шестигранника, посчитать матрицу Якоби и т.д.
//-------------------------------------------------------------------------------
T_Brick::T_Brick(double *x_coords, double *y_coords, double *z_coords, long type_of_hex)
{
	for (long i=0; i<8; i++)
	{
		this->x[i] = x_coords[i];
		this->y[i] = y_coords[i];
		this->z[i] = z_coords[i];
	}

	this->xk  = this->x[0];
	this->xk1 = this->x[7];
	this->yk  = this->y[0];
	this->yk1 = this->y[7];
	this->zk  = this->z[0];
	this->zk1 = this->z[7];

	this->hx = this->xk1 - this->xk;
	this->hy = this->yk1 - this->yk; 
	this->hz = this->zk1 - this->zk;

	this->type_of_hex = type_of_hex;
}
//------------------------------------------------------------------------
T_Brick::T_Brick(long num, long (*nver)[14], double (*xyz)[3])
{
	this->num = num;
	this->nver = nver;
	this->xyz = xyz;

	for (long i=0; i<8; i++)
	{
		this->x[i] = xyz[nver[num][i]][0];
		this->y[i] = xyz[nver[num][i]][1];
		this->z[i] = xyz[nver[num][i]][2];
	}

	this->type_of_hex = nver[num][13];

	if (this->type_of_hex <= 30)
	{
		this->xk  = this->x[0];
		this->xk1 = this->x[7];
		this->yk  = this->y[0];
		this->yk1 = this->y[7];
		this->zk  = this->z[0];
		this->zk1 = this->z[7];

		this->hx = this->xk1 - this->xk;
		this->hy = this->yk1 - this->yk; 
		this->hz = this->zk1 - this->zk;
	}
}
//-------------------------------------------------------------------------------
//------------    конструктор для МТЗ   -----------------------------------------
//-------------------------------------------------------------------------------
T_Brick::T_Brick(long num, long (*nver)[14], long (*ed)[25], long (*edges)[2], double (*xyz)[3], long *nvkat,
	double *sigma3d, double *sigma0, double *mu3d,  double omega,
	long alpha, long n_1d, double *z_1d, double *sin_1d, double *cos_1d)
{
	this->num = num;
	this->nver = nver;
	this->ed = ed;
	this->xyz = xyz;
	this->edges = edges;

	for (long i=0; i<8; i++)
	{
		this->x[i] = xyz[nver[num][i]][0];
		this->y[i] = xyz[nver[num][i]][1];
		this->z[i] = xyz[nver[num][i]][2];
	}

	this->xk  = this->x[0];
	this->xk1 = this->x[7];
	this->yk  = this->y[0];
	this->yk1 = this->y[7];
	this->zk  = this->z[0];
	this->zk1 = this->z[7];

	this->hx = this->xk1 - this->xk;
	this->hy = this->yk1 - this->yk; 
	this->hz = this->zk1 - this->zk;

	this->nvkat = nvkat;
	this->n_mat = nvkat[num];
	this->sigma = sigma3d[nvkat[num]];
	this->sigma0 = sigma0[nvkat[num]];
	this->mu = mu3d[nvkat[num]];
	this->omega = omega;

	this->mu0 = this->mu;
	this->dpr0 = this->dpr = 0;

	this->type_of_hex = nver[num][13];

	this->alpha = alpha;
	this->n_1d =n_1d;
	this->z_1d = z_1d;
	this->sin_1d = sin_1d;
	this->cos_1d = cos_1d;
}
//-------------------------------------------------------------------------------
//----  конструктор для нестационарной задачи
//-------------------------------------------------------------------------------
T_Brick::T_Brick(long num, long (*nver)[14], long (*ed)[25], long (*edges)[2], double (*xyz)[3], long *nvkat,
		double *sigma_table,  double *sigma0_table, double *mu_table, 
		double *En)
{
	this->num = num;
	this->nver = nver;
	this->ed = ed;
	this->xyz = xyz;
	this->edges = edges;

	for (long i=0; i<8; i++)
	{
		this->x[i] = xyz[nver[num][i]][0];
		this->y[i] = xyz[nver[num][i]][1];
		this->z[i] = xyz[nver[num][i]][2];
	}

	this->xk  = this->x[0];
	this->xk1 = this->x[7];
	this->yk  = this->y[0];
	this->yk1 = this->y[7];
	this->zk  = this->z[0];
	this->zk1 = this->z[7];

	this->hx = this->xk1 - this->xk;
	this->hy = this->yk1 - this->yk; 
	this->hz = this->zk1 - this->zk;

	this->nvkat = nvkat;
	this->sigma = sigma_table[nvkat[num]];
	this->sigma0 = sigma0_table[nvkat[num]];
	this->mu = mu_table[nvkat[num]];

	this->mu0 = this->mu;
	this->dpr0 = this->dpr = 0;

	this->En = En;

	this->type_of_hex = nver[num][13];

	this->n_mat = nvkat[num];
}
//-------------------------------------------------------------------------------
//------  Деструктор
//-------------------------------------------------------------------------------
T_Brick::~T_Brick()
{
}
void T_Brick::Calc_ss0()
{
	this->s_s0 = this->sigmaTensor;
	for(int ii=0;ii<3;ii++)
	{
		for(int jj=0;jj<3;jj++)
		{
			this->s_s0.val[ii][jj] -= this->sigmaTensor0.val[ii][jj];
		}
	}
}
//-------------------------------------------------------------------------------
//------- Вычислить локальную матрицу жёсткости на параллелепипеде --------------
//-------------------------------------------------------------------------------
void T_Brick::Compute_Local_Matrix_B()
{
	double mult_x, mult_y, mult_z;
	double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15;

	mult_x = 2.0*hx/(hy*hz);
	mult_y = 2.0*hy/(hx*hz);
	mult_z = 2.0*hz/(hx*hy);

	t1 = mult_y+mult_z;
	t2 = mult_y/3.0;
	t3 = 2.0/3.0*mult_z;
	t4 = t2-t3;
	t5 = 2.0/3.0*mult_y;
	t6 = mult_z/3.0;
	t7 = -t5+t6;
	t8 = mult_x+mult_z;
	t9 = 2.0/3.0*mult_x;
	t10 = -t9+t6;
	t11 = mult_x/3.0;
	t12 = t11-t3;
	t13 = mult_x+mult_y;
	t14 = t11-t5;
	t15 = -t9+t2;
	b[0][0] = 2.0/3.0*t1;
	b[0][1] = t4;
	b[0][2] = t7;
	b[0][3] = -t1/3.0;
	b[0][4] = -t3;
	b[0][5] = -t6;
	b[0][6] = t3;
	b[0][7] = t6;
	b[0][8] = -t5;
	b[0][9] = t5;
	b[0][10] = -t2;
	b[0][11] = t2;
	b[1][0] = t4;
	b[1][1] = 2.0/3.0*t1;
	b[1][2] = -t1/3.0;
	b[1][3] = t7;
	b[1][4] = t3;
	b[1][5] = t6;
	b[1][6] = -t3;
	b[1][7] = -t6;
	b[1][8] = -t2;
	b[1][9] = t2;
	b[1][10] = -t5;
	b[1][11] = t5;
	b[2][0] = t7;
	b[2][1] = -t1/3.0;
	b[2][2] = 2.0/3.0*t1;
	b[2][3] = t4;
	b[2][4] = -t6;
	b[2][5] = -t3;
	b[2][6] = t6;
	b[2][7] = t3;
	b[2][8] = t5;
	b[2][9] = -t5;
	b[2][10] = t2;
	b[2][11] = -t2;
	b[3][0] = -t1/3.0;
	b[3][1] = t7;
	b[3][2] = t4;
	b[3][3] = 2.0/3.0*t1;
	b[3][4] = t6;
	b[3][5] = t3;
	b[3][6] = -t6;
	b[3][7] = -t3;
	b[3][8] = t2;
	b[3][9] = -t2;
	b[3][10] = t5;
	b[3][11] = -t5;
	b[4][0] = -t3;
	b[4][1] = t3;
	b[4][2] = -t6;
	b[4][3] = t6;
	b[4][4] = 2.0/3.0*t8;
	b[4][5] = t10;
	b[4][6] = t12;
	b[4][7] = -t8/3.0;
	b[4][8] = -t9;
	b[4][9] = -t11;
	b[4][10] = t9;
	b[4][11] = t11;
	b[5][0] = -t6;
	b[5][1] = t6;
	b[5][2] = -t3;
	b[5][3] = t3;
	b[5][4] = t10;
	b[5][5] = 2.0/3.0*t8;
	b[5][6] = -t8/3.0;
	b[5][7] = t12;
	b[5][8] = t9;
	b[5][9] = t11;
	b[5][10] = -t9;
	b[5][11] = -t11;
	b[6][0] = t3;
	b[6][1] = -t3;
	b[6][2] = t6;
	b[6][3] = -t6;
	b[6][4] = t12;
	b[6][5] = -t8/3.0;
	b[6][6] = 2.0/3.0*t8;
	b[6][7] = t10;
	b[6][8] = -t11;
	b[6][9] = -t9;
	b[6][10] = t11;
	b[6][11] = t9;
	b[7][0] = t6;
	b[7][1] = -t6;
	b[7][2] = t3;
	b[7][3] = -t3;
	b[7][4] = -t8/3.0;
	b[7][5] = t12;
	b[7][6] = t10;
	b[7][7] = 2.0/3.0*t8;
	b[7][8] = t11;
	b[7][9] = t9;
	b[7][10] = -t11;
	b[7][11] = -t9;
	b[8][0] = -t5;
	b[8][1] = -t2;
	b[8][2] = t5;
	b[8][3] = t2;
	b[8][4] = -t9;
	b[8][5] = t9;
	b[8][6] = -t11;
	b[8][7] = t11;
	b[8][8] = 2.0/3.0*t13;
	b[8][9] = t14;
	b[8][10] = t15;
	b[8][11] = -t13/3.0;
	b[9][0] = t5;
	b[9][1] = t2;
	b[9][2] = -t5;
	b[9][3] = -t2;
	b[9][4] = -t11;
	b[9][5] = t11;
	b[9][6] = -t9;
	b[9][7] = t9;
	b[9][8] = t14;
	b[9][9] = 2.0/3.0*t13;
	b[9][10] = -t13/3.0;
	b[9][11] = t15;
	b[10][0] = -t2;
	b[10][1] = -t5;
	b[10][2] = t2;
	b[10][3] = t5;
	b[10][4] = t9;
	b[10][5] = -t9;
	b[10][6] = t11;
	b[10][7] = -t11;
	b[10][8] = t15;
	b[10][9] = -t13/3.0;
	b[10][10] = 2.0/3.0*t13;
	b[10][11] = t14;
	b[11][0] = t2;
	b[11][1] = t5;
	b[11][2] = -t2;
	b[11][3] = -t5;
	b[11][4] = t11;
	b[11][5] = -t11;
	b[11][6] = t9;
	b[11][7] = -t9;
	b[11][8] = -t13/3.0;
	b[11][9] = t15;
	b[11][10] = t14;
	b[11][11] = 2.0/3.0*t13;
}
//-------------------------------------------------------------------------------
void T_Brick::Compute_Local_Matrix_C()
{
	{
		double k11=sigmaTensor.val[0][0]*hy*hz/9./hx;
		double k22=sigmaTensor.val[1][1]*hx*hz/9./hy;
		double k33=sigmaTensor.val[2][2]*hx*hy/9./hz;
		double k21=sigmaTensor.val[1][0]*hz/6.;
		double k31=sigmaTensor.val[2][0]*hy/6.;
		double k32=sigmaTensor.val[2][1]*hx/6.;

		double k114=k11*4;
		double k112=k11*2;
		double k212=k21*2;
		double k312=k31*2;
		double k224=k22*4;
		double k222=k22*2;
		double k322=k32*2;
		double k334=k33*4;
		double k332=k33*2;

		double C[12][12]=
		{
			{k114,k112,k112,k11,	k212,k212,k21,k21,	k312,k312,k31,k31},
			{k112,k114,k11,k112,	k212,k212,k21,k21,	k31,k31,k312,k312},
			{k112,k11,k114,k112,	k21,k21,k212,k212,	k312,k312,k31,k31},
			{k11,k112,k112,k114,	k21,k21,k212,k212,	k31,k31,k312,k312},

			{k212,k212,k21,k21,		k224,k222,k222,k22,	k322,k32,k322,k32},
			{k212,k212,k21,k21,		k222,k224,k22,k222,	k32,k322,k32,k322},
			{k21,k21,k212,k212,		k222,k22,k224,k222,	k322,k32,k322,k32},
			{k21,k21,k212,k212,		k22,k222,k222,k224,	k32,k322,k32,k322},

			{k312,k31,k312,k31,		k322,k32,k322,k32,	k334,k332,k332,k33},
			{k312,k31,k312,k31,		k32,k322,k32,k322,	k332,k334,k33,k332},
			{k31,k312,k31,k312,		k322,k32,k322,k32,	k332,k33,k334,k332},
			{k31,k312,k31,k312,		k32,k322,k32,k322,	k33,k332,k332,k334}
		};

		for (int i=0;i<12;i++)
		{
			for (int j=0;j<12;j++)
			{
				cSigma[i][j]=C[i][j];
			}
		}
	}
	{
		double k11=(s_s0.val[0][0])*hy*hz/9./hx;
		double k22=(s_s0.val[1][1])*hx*hz/9./hy;
		double k33=(s_s0.val[2][2])*hx*hy/9./hz;
		double k21=(s_s0.val[1][0])*hz/6.;
		double k31=(s_s0.val[2][0])*hy/6.;
		double k32=(s_s0.val[2][1])*hx/6.;

		double k114=k11*4;
		double k112=k11*2;
		double k212=k21*2;
		double k312=k31*2;
		double k224=k22*4;
		double k222=k22*2;
		double k322=k32*2;
		double k334=k33*4;
		double k332=k33*2;

		double C[12][12]=
		{
			{k114,k112,k112,k11,	k212,k212,k21,k21,	k312,k312,k31,k31},
			{k112,k114,k11,k112,	k212,k212,k21,k21,	k31,k31,k312,k312},
			{k112,k11,k114,k112,	k21,k21,k212,k212,	k312,k312,k31,k31},
			{k11,k112,k112,k114,	k21,k21,k212,k212,	k31,k31,k312,k312},

			{k212,k212,k21,k21,		k224,k222,k222,k22,	k322,k32,k322,k32},
			{k212,k212,k21,k21,		k222,k224,k22,k222,	k32,k322,k32,k322},
			{k21,k21,k212,k212,		k222,k22,k224,k222,	k322,k32,k322,k32},
			{k21,k21,k212,k212,		k22,k222,k222,k224,	k32,k322,k32,k322},

			{k312,k31,k312,k31,		k322,k32,k322,k32,	k334,k332,k332,k33},
			{k312,k31,k312,k31,		k32,k322,k32,k322,	k332,k334,k33,k332},
			{k31,k312,k31,k312,		k322,k32,k322,k32,	k332,k33,k334,k332},
			{k31,k312,k31,k312,		k32,k322,k32,k322,	k33,k332,k332,k334}
		};

		for (int i=0;i<12;i++)
		{
			for (int j=0;j<12;j++)
			{
				cSigma0[i][j]=C[i][j];
			}
		}
	}
	{
		double k11=hy*hz/9./hx;
		double k22=hx*hz/9./hy;
		double k33=hx*hy/9./hz;
		double k21=0;
		double k31=0;
		double k32=0;

		double k114=k11*4;
		double k112=k11*2;
		double k212=k21*2;
		double k312=k31*2;
		double k224=k22*4;
		double k222=k22*2;
		double k322=k32*2;
		double k334=k33*4;
		double k332=k33*2;

		double C[12][12]=
		{
			{k114,k112,k112,k11,	k212,k212,k21,k21,	k312,k312,k31,k31},
			{k112,k114,k11,k112,	k212,k212,k21,k21,	k31,k31,k312,k312},
			{k112,k11,k114,k112,	k21,k21,k212,k212,	k312,k312,k31,k31},
			{k11,k112,k112,k114,	k21,k21,k212,k212,	k31,k31,k312,k312},

			{k212,k212,k21,k21,		k224,k222,k222,k22,	k322,k32,k322,k32},
			{k212,k212,k21,k21,		k222,k224,k22,k222,	k32,k322,k32,k322},
			{k21,k21,k212,k212,		k222,k22,k224,k222,	k322,k32,k322,k32},
			{k21,k21,k212,k212,		k22,k222,k222,k224,	k32,k322,k32,k322},

			{k312,k31,k312,k31,		k322,k32,k322,k32,	k334,k332,k332,k33},
			{k312,k31,k312,k31,		k32,k322,k32,k322,	k332,k334,k33,k332},
			{k31,k312,k31,k312,		k322,k32,k322,k32,	k332,k33,k334,k332},
			{k31,k312,k31,k312,		k32,k322,k32,k322,	k33,k332,k332,k334}
		};

		for (int i=0;i<12;i++)
		{
			for (int j=0;j<12;j++)
			{
				c0[i][j]=C[i][j];
			}
		}
	}
}
//---------------------------------------------------------------------------------
//-------   Вычислить локальную матрицу и вектор для нестационарной задачии -------
//---------------------------------------------------------------------------------
void T_Brick::Compute_Local_Matrix_And_Vector(const long what_compute)
{
	long i,j;

	switch(what_compute)
	{
	case 1: // вычислить матрицу массы и вектор
		Calc_local_matrix_c_for_hexahedron();

		for(i=0; i<12; i++)
		{
			for(j=0; j<12; j++)
				a[i][j] = c[i][j]*sigma;
			g[i] = 0.0;
		}
		break;

	case 2: // если требуется только прибавить матрицу жёсткости
		Calc_local_matrix_b_for_hexahedron();

		for(i=0; i<12; i++)
		{
			for(j=0; j<12; j++)
				a[i][j] = b[i][j]/mu;
			g[i] = 0.0;
		}
		break;

//  [24/3/2008 Domnikov]
	case 3: // собрать матрицу массы с коэффициентом диэлектрической проницаемости
		Calc_local_matrix_c_for_hexahedron();

		for(i=0; i<12; i++)
		{
			for(j=0; j<12; j++)
				a[i][j] = c[i][j]*dpr;
			g[i] = 0.0;
		}

		break;
	}
}
//-------------------------------------------------------------------------------
//-----   Вычислить локальный вектор правой части для нестационарной задачи   ---
//-----   через En (нормальное поле, заданное на рёбрах)   ----------------------
//-------------------------------------------------------------------------------
void T_Brick::Compute_Local_Vector_For_Anomal_Problem()
{
	long i;
	double temp;
	long edge;

	temp = sigma - sigma0;

	for(i=0; i<12; i++)
	{
		edge = ed[num][i];
		f_re[i] = En[edge]*temp;
	}

	Mult_Plot((double*)c, f_re, g, 12);
}
//------------------------------------------------------------------------
// вычисление вектора правой части через нормальное поле, в случае, когда все коэф-ты разрывны
//------------------------------------------------------------------------
void T_Brick::ComputeLocalVectorMuEpsSigma(double *An, double *d2An)
{
	int i;
	double t_sig, t_eps, t_mu;
	long edge;
	double g_sig[12], g_mu[12];
	double f_sig[12], f_mu[12];

	t_sig = sigma - sigma0;
	t_eps = dpr0 - dpr;
	t_mu = 1/mu0 - 1/mu;

	for(i=0; i<12; i++)
	{
		edge = ed[num][i];
		f_sig[i] = En[edge]*t_sig + d2An[edge]*t_eps;
		f_mu[i]  = An[edge]*t_mu;
	}

	Mult_Plot((double*)b, f_mu, g_mu, 12);
	Mult_Plot((double*)c, f_sig, g_sig, 12);

	for(i=0; i<12; i++)
		g[i] = g_mu[i] + g_sig[i];
}
//-------------------------------------------------------------------------------
//---- Вычислить матрицу Якоби в точке Гаусса с номером n_of_point  -------------
//-------------------------------------------------------------------------------
void T_Brick::Calc_J(int n_of_point)
{
	long i, j;

	if (type_of_hex > 30) //  шестигранник
	{
		for(i=0; i<3; i++)
			for(j=0; j<3; j++)
				J[i][j] = 0.0;

		// эл-ты матрицы Якоби
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

		// вычисляем Якобиан (определитель)
		det_J = J[0][0]*J[1][1]*J[2][2] - J[0][0]*J[1][2]*J[2][1] + J[1][0]*J[2][1]*J[0][2]
		- J[1][0]*J[0][1]*J[2][2] + J[2][0]*J[0][1]*J[1][2] - J[2][0]*J[1][1]*J[0][2];

		// модуль Якобиана
		det_J_abs = fabs(det_J);

		// матрица, обратная к матрице Якоби (и транспонированная)
		J_1_T[0][0] = J_1[0][0] = (J[1][1]*J[2][2]-J[2][1]*J[1][2])/det_J;
		J_1_T[1][0] = J_1[0][1] = (J[2][1]*J[0][2]-J[0][1]*J[2][2])/det_J;
		J_1_T[2][0] = J_1[0][2] = (J[0][1]*J[1][2]-J[1][1]*J[0][2])/det_J;
		J_1_T[0][1] = J_1[1][0] = (-J[1][0]*J[2][2]+J[2][0]*J[1][2])/det_J;
		J_1_T[1][1] = J_1[1][1] = (J[0][0]*J[2][2]-J[2][0]*J[0][2])/det_J;
		J_1_T[2][1] = J_1[1][2] = (-J[0][0]*J[1][2]+J[1][0]*J[0][2])/det_J;
		J_1_T[0][2] = J_1[2][0] = (J[1][0]*J[2][1]-J[2][0]*J[1][1])/det_J;
		J_1_T[1][2] = J_1[2][1] = (-J[0][0]*J[2][1]+J[2][0]*J[0][1])/det_J;
		J_1_T[2][2] = J_1[2][2] = (J[0][0]*J[1][1]-J[1][0]*J[0][1])/det_J;
	} 
	else // параллелепипед
	{
		Calc_J_in_parallelepiped();		
	}
}
//-------------------------------------------------------------------------------
//---- Вычислить матрицу Якоби в произвольной точке внутри шестигранника --------
//-------------------------------------------------------------------------------
void T_Brick::Calc_J(double xi, double eta, double zeta)
{
	long i, j;
	double dphi_x, dphi_y, dphi_z;

	if (type_of_hex > 30) // шестигранник
	{
		for(i=0; i<3; i++)
			for(j=0; j<3; j++)
				J[i][j] = 0.0;

		// эл-ты матрицы Якоби
		for(i=0; i<8; i++)
		{
			dphi_x = dPhi_node(i, 0, xi, eta, zeta);
			dphi_y = dPhi_node(i, 1, xi, eta, zeta);
			dphi_z = dPhi_node(i, 2, xi, eta, zeta);

			J[0][0] += x[i]*dphi_x; //  d_xi[i];
			J[0][1] += x[i]*dphi_y; //  d_eta[i];
			J[0][2] += x[i]*dphi_z; //  d_zeta[i];

			J[1][0] += y[i]*dphi_x; //  d_xi[i];
			J[1][1] += y[i]*dphi_y; //  d_eta[i];
			J[1][2] += y[i]*dphi_z; //  d_zeta[i];

			J[2][0] += z[i]*dphi_x; //  d_xi[i];
			J[2][1] += z[i]*dphi_y; //  d_eta[i];
			J[2][2] += z[i]*dphi_z; //  d_zeta[i];
		}

		// вычисляем Якобиан (определитель)
		det_J = J[0][0]*J[1][1]*J[2][2] - J[0][0]*J[1][2]*J[2][1] + J[1][0]*J[2][1]*J[0][2]
		- J[1][0]*J[0][1]*J[2][2] + J[2][0]*J[0][1]*J[1][2] - J[2][0]*J[1][1]*J[0][2];

		// модуль Якобиана
		det_J_abs = fabs(det_J);

		// матрица, обратная к матрице Якоби (и транспонированная)
		J_1_T[0][0] = J_1[0][0] = (J[1][1]*J[2][2]-J[2][1]*J[1][2])/det_J;
		J_1_T[1][0] = J_1[0][1] = (J[2][1]*J[0][2]-J[0][1]*J[2][2])/det_J;
		J_1_T[2][0] = J_1[0][2] = (J[0][1]*J[1][2]-J[1][1]*J[0][2])/det_J;
		J_1_T[0][1] = J_1[1][0] = (-J[1][0]*J[2][2]+J[2][0]*J[1][2])/det_J;
		J_1_T[1][1] = J_1[1][1] = (J[0][0]*J[2][2]-J[2][0]*J[0][2])/det_J;
		J_1_T[2][1] = J_1[1][2] = (-J[0][0]*J[1][2]+J[1][0]*J[0][2])/det_J;
		J_1_T[0][2] = J_1[2][0] = (J[1][0]*J[2][1]-J[2][0]*J[1][1])/det_J;
		J_1_T[1][2] = J_1[2][1] = (-J[0][0]*J[2][1]+J[2][0]*J[0][1])/det_J;
		J_1_T[2][2] = J_1[2][2] = (J[0][0]*J[1][1]-J[1][0]*J[0][1])/det_J;

	} 
	else // параллелепипед
	{
		Calc_J_in_parallelepiped();
	}
}

void T_Brick::Calc_V_Node_2(double *q,double xi, double eta, double zeta)
{
	long i;
	double dphi_x, dphi_y, dphi_z;

	for(i=0; i<3; i++)V[i] = 0.0;

	for(i=0; i<8; i++)
	{
		dphi_x = dPhi_node_2(i, 0, xi, eta, zeta);
		dphi_y = dPhi_node_2(i, 1, xi, eta, zeta);
		dphi_z = dPhi_node_2(i, 2, xi, eta, zeta);
		
		V[0] += q[i]*dphi_x; //  d_xi[i];
		V[1] += q[i]*dphi_y; //  d_eta[i];
		V[2] += q[i]*dphi_z; //  d_zeta[i];
	}
}

void T_Brick::Calc_J_Node_2(double xi, double eta, double zeta)
{
	long i, j;
	double dphi_x, dphi_y, dphi_z;

	if (type_of_hex > 30) // шестигранник
	{
		for(i=0; i<3; i++)
			for(j=0; j<3; j++)
				J[i][j] = 0.0;

		// эл-ты матрицы Якоби
		for(i=0; i<8; i++)
		{
			dphi_x = dPhi_node_2(i, 0, xi, eta, zeta);
			dphi_y = dPhi_node_2(i, 1, xi, eta, zeta);
			dphi_z = dPhi_node_2(i, 2, xi, eta, zeta);
			
			J[0][0] += x[i]*dphi_x; //  d_xi[i];
			J[0][1] += x[i]*dphi_y; //  d_eta[i];
			J[0][2] += x[i]*dphi_z; //  d_zeta[i];
			
			J[1][0] += y[i]*dphi_x; //  d_xi[i];
			J[1][1] += y[i]*dphi_y; //  d_eta[i];
			J[1][2] += y[i]*dphi_z; //  d_zeta[i];
			
			J[2][0] += z[i]*dphi_x; //  d_xi[i];
			J[2][1] += z[i]*dphi_y; //  d_eta[i];
			J[2][2] += z[i]*dphi_z; //  d_zeta[i];
		}
		// вычисляем Якобиан (определитель)
		det_J = J[0][0]*J[1][1]*J[2][2] - J[0][0]*J[1][2]*J[2][1] + J[1][0]*J[2][1]*J[0][2]
		- J[1][0]*J[0][1]*J[2][2] + J[2][0]*J[0][1]*J[1][2] - J[2][0]*J[1][1]*J[0][2];
		// модуль Якобиана
		det_J_abs = fabs(det_J);
		// матрица, обратная к матрице Якоби (и транспонированная)
		J_1_T[0][0] = J_1[0][0] = (J[1][1]*J[2][2]-J[2][1]*J[1][2])/det_J;
		J_1_T[1][0] = J_1[0][1] = (J[2][1]*J[0][2]-J[0][1]*J[2][2])/det_J;
		J_1_T[2][0] = J_1[0][2] = (J[0][1]*J[1][2]-J[1][1]*J[0][2])/det_J;
		J_1_T[0][1] = J_1[1][0] = (-J[1][0]*J[2][2]+J[2][0]*J[1][2])/det_J;
		J_1_T[1][1] = J_1[1][1] = (J[0][0]*J[2][2]-J[2][0]*J[0][2])/det_J;
		J_1_T[2][1] = J_1[1][2] = (-J[0][0]*J[1][2]+J[1][0]*J[0][2])/det_J;
		J_1_T[0][2] = J_1[2][0] = (J[1][0]*J[2][1]-J[2][0]*J[1][1])/det_J;
		J_1_T[1][2] = J_1[2][1] = (-J[0][0]*J[2][1]+J[2][0]*J[0][1])/det_J;
		J_1_T[2][2] = J_1[2][2] = (J[0][0]*J[1][1]-J[1][0]*J[0][1])/det_J;

	} 
	else // параллелепипед
	{
		// эл-ты матрицы Якоби
		J[0][0] = hx;
		J[1][1] = hy;
		J[2][2] = hz;
		J[0][1] = J[0][2] = J[1][0] = J[1][2] = J[2][0] = J[2][1] = 0.0;
		// вычисляем Якобиан (определитель)
		det_J = hx*hy*hz;
		// модуль Якобиана
		det_J_abs = fabs(det_J);
		// матрица, обратная к матрице Якоби (и транспонированная)
		J_1_T[0][0] = J_1[0][0] = 1.0/hx;
		J_1_T[1][1] = J_1[1][1] = 1.0/hy;
		J_1_T[2][2] = J_1[2][2] = 1.0/hz;
		J_1_T[0][1] = J_1_T[0][2] = J_1_T[1][0] = J_1_T[1][2] = J_1_T[2][0] = J_1_T[2][1] =
		J_1[0][1] = J_1[0][2] = J_1[1][0] = J_1[1][2] = J_1[2][0] = J_1[2][1] =	0.0;
	}
}

//-----------------------------------------------------------------------------
//----- Вычислить матрицу Якоби в точках Гаусса, расположенных на верхней -----
//----- грани шестигранника (это нужно для выдачи) ----------------------------
//-----------------------------------------------------------------------------
void T_Brick::Calc_J_on_face(int n_of_point)
{
	long i, j;

	if (type_of_hex > 30) // шестигранник
	{
		for(i=0; i<3; i++)
			for(j=0; j<3; j++)
				J[i][j] = 0.0;

		// эл-ты матрицы Якоби
		for(i=0; i<8; i++)
		{
			J[0][0] += x[i]*gauss_3_d_phi_face[n_of_point][i][0]; //  d_xi[i];
			J[0][1] += x[i]*gauss_3_d_phi_face[n_of_point][i][1]; //  d_eta[i];
			J[0][2] += x[i]*gauss_3_d_phi_face[n_of_point][i][2]; //  d_zeta[i];

			J[1][0] += y[i]*gauss_3_d_phi_face[n_of_point][i][0]; //  d_xi[i];
			J[1][1] += y[i]*gauss_3_d_phi_face[n_of_point][i][1]; //  d_eta[i];
			J[1][2] += y[i]*gauss_3_d_phi_face[n_of_point][i][2]; //  d_zeta[i];

			J[2][0] += z[i]*gauss_3_d_phi_face[n_of_point][i][0]; //  d_xi[i];
			J[2][1] += z[i]*gauss_3_d_phi_face[n_of_point][i][1]; //  d_eta[i];
			J[2][2] += z[i]*gauss_3_d_phi_face[n_of_point][i][2]; //  d_zeta[i];
		}

		// вычисляем Якобиан (определитель)
		det_J = J[0][0]*J[1][1]*J[2][2] - J[0][0]*J[1][2]*J[2][1] + J[1][0]*J[2][1]*J[0][2]
		- J[1][0]*J[0][1]*J[2][2] + J[2][0]*J[0][1]*J[1][2] - J[2][0]*J[1][1]*J[0][2];

		// модуль Якобиана
		det_J_abs = fabs(det_J);

		// матрица, обратная к матрице Якоби (и транспонированная)
		J_1_T[0][0] = J_1[0][0] = (J[1][1]*J[2][2]-J[2][1]*J[1][2])/det_J;
		J_1_T[1][0] = J_1[0][1] = (J[2][1]*J[0][2]-J[0][1]*J[2][2])/det_J;
		J_1_T[2][0] = J_1[0][2] = (J[0][1]*J[1][2]-J[1][1]*J[0][2])/det_J;
		J_1_T[0][1] = J_1[1][0] = (-J[1][0]*J[2][2]+J[2][0]*J[1][2])/det_J;
		J_1_T[1][1] = J_1[1][1] = (J[0][0]*J[2][2]-J[2][0]*J[0][2])/det_J;
		J_1_T[2][1] = J_1[1][2] = (-J[0][0]*J[1][2]+J[1][0]*J[0][2])/det_J;
		J_1_T[0][2] = J_1[2][0] = (J[1][0]*J[2][1]-J[2][0]*J[1][1])/det_J;
		J_1_T[1][2] = J_1[2][1] = (-J[0][0]*J[2][1]+J[2][0]*J[0][1])/det_J;
		J_1_T[2][2] = J_1[2][2] = (J[0][0]*J[1][1]-J[1][0]*J[0][1])/det_J;
	} 
	else // параллелепипед
	{
		Calc_J_in_parallelepiped();
	}
}
//---------------------------------------------------------------------------
//------ Одномерные шаблонные базисные функции на [-1, 1]  ------------------
//---------------------------------------------------------------------------
double T_Brick::l0(double x)
{
	return (1.0 - x)*0.5;
}
//--------------------------------
double T_Brick::l1(double x)
{
	return (x + 1.0)*0.5;
}
//-----------------------------------------------------------------------------------------
//---- Узловые базисные функции на параллелепипеде [-1,1]^3 (нужны для преобразования) ----
//-----------------------------------------------------------------------------------------
double T_Brick::Phi_node(long i, double x, double y, double z)
{
	switch(i)
	{
	case 0:
		return l0(x)*l0(y)*l0(z);
	case 1:
		return l1(x)*l0(y)*l0(z);
	case 2:
		return l0(x)*l1(y)*l0(z);
	case 3:
		return l1(x)*l1(y)*l0(z);
	case 4:
		return l0(x)*l0(y)*l1(z);
	case 5:
		return l1(x)*l0(y)*l1(z);
	case 6:
		return l0(x)*l1(y)*l1(z);
	default:
		return l1(x)*l1(y)*l1(z);
	}
}

double T_Brick::dPhi_node_2(long i, long j, double xi, double eta, double zeta)
{
	switch(j)
	{
	case 0:
		switch(i)
		{
		case 0: return -(1.0-eta)*(1.0-zeta);
		case 1: return  (1.0-eta)*(1.0-zeta);
		case 2: return -(eta)*(1.0-zeta);
		case 3: return  (eta)*(1.0-zeta);
		case 4: return -(1.0-eta)*(zeta);
		case 5: return  (1.0-eta)*(zeta);
		case 6: return -(eta)*(zeta);
		default: return (eta)*(zeta);
		}
	case 1:
		switch(i)
		{
		case 0: return -(1.0-xi)*(1.0-zeta);
		case 1: return -(xi)*(1.0-zeta);
		case 2: return (1.0-xi)*(1.0-zeta);
		case 3: return (xi)*(1.0-zeta);
		case 4: return -(1.0-xi)*(zeta);
		case 5: return -(xi)*(zeta);
		case 6: return  (1.0-xi)*(zeta);
		default: return (xi)*(zeta);
		}
	default:
		switch(i)
		{
		case 0: return -(1.0-xi)*(1.0-eta);
		case 1: return -(xi)*(1.0-eta);
		case 2: return -(1.0-xi)*(eta);
		case 3: return -(xi)*(eta);
		case 4: return  (1.0-xi)*(1.0-eta);
		case 5: return  (xi)*(1.0-eta);
		case 6: return  (1.0-xi)*(eta);
		default: return (xi)*(eta);
		}
	}
}
//-------------------------------------------------------------------------------
//--- Производные от узловых базисных функций на параллелепипеде [-1, 1]^3 ------
//--- (нужны для преобразования) ------------------------------------------------
//-------------------------------------------------------------------------------
double T_Brick::dPhi_node(long i, long j, double xi, double eta, double zeta)
{
	switch(j)
	{
	case 0:
		switch(i)
		{
		case 0: return -(-1.0+eta)*(-1.0+zeta)/8.0;
		case 1: return (-1.0+eta)*(-1.0+zeta)/8.0;
		case 2: return (eta+1.0)*(-1.0+zeta)/8.0;
		case 3: return -(eta+1.0)*(-1.0+zeta)/8.0;
		case 4: return (-1.0+eta)*(zeta+1.0)/8.0;
		case 5: return -(-1.0+eta)*(zeta+1.0)/8.0;
		case 6: return -(eta+1.0)*(zeta+1.0)/8.0;
		default: return (eta+1.0)*(zeta+1.0)/8.0;
		}
	case 1:
		switch(i)
		{
		case 0: return -(-1.0+xi)*(-1.0+zeta)/8.0;      
		case 1: return (xi+1.0)*(-1.0+zeta)/8.0;
		case 2: return (-1.0+xi)*(-1.0+zeta)/8.0;
		case 3: return -(xi+1.0)*(-1.0+zeta)/8.0;
		case 4: return (-1.0+xi)*(zeta+1.0)/8.0;
		case 5: return -(xi+1.0)*(zeta+1.0)/8.0;
		case 6: return -(-1.0+xi)*(zeta+1.0)/8.0;
		default: return (xi+1.0)*(zeta+1.0)/8.0;
		}
	default:
		switch(i)
		{
		case 0: return -(-1.0+xi)*(-1.0+eta)/8.0;      
		case 1: return (xi+1.0)*(-1.0+eta)/8.0;
		case 2: return (-1.0+xi)*(eta+1.0)/8.0;
		case 3: return -(xi+1.0)*(eta+1.0)/8.0;
		case 4: return (-1.0+xi)*(-1.0+eta)/8.0;
		case 5: return -(xi+1.0)*(-1.0+eta)/8.0;
		case 6: return -(-1.0+xi)*(eta+1.0)/8.0;
		default: return (xi+1.0)*(eta+1.0)/8.0;
		}
	}
}
//-----------------------------------------------------------------------------
//---- Шаблонные базисные функции на векторном параллелепипеде [-1,1]^3 -------
//-----------------------------------------------------------------------------
void T_Brick::Basis_func_on_reference_vec_par(long n, double *in, double *out)
{
	double x, y, z;

	x = in[0];
	y = in[1];
	z = in[2];

	switch(n)
	{
	case 0:
		out[0] = l0(y)*l0(z);
		out[1] = 0.0;
		out[2] = 0.0;
		break;
	case 1:
		out[0] = l1(y)*l0(z); 
		out[1] = 0.0;
		out[2] = 0.0;
		break;
	case 2:
		out[0] = l0(y)*l1(z);
		out[1] = 0.0;
		out[2] = 0.0;
		break;
	case 3:
		out[0] = l1(y)*l1(z);
		out[1] = 0.0;
		out[2] = 0.0;
		break;
	case 4:
		out[0] = 0.0;
		out[1] = l0(z)*l0(x);
		out[2] = 0.0;
		break;
	case 5:
		out[0] = 0.0;
		out[1] = l1(z)*l0(x);
		out[2] = 0.0;
		break;
	case 6:
		out[0] = 0.0;
		out[1] = l0(z)*l1(x);
		out[2] = 0.0;
		break;
	case 7:
		out[0] = 0.0;
		out[1] = l1(z)*l1(x);
		out[2] = 0.0;
		break;
	case 8:
		out[0] = 0.0;
		out[1] = 0.0;
		out[2] = l0(x)*l0(y);
		break;
	case 9:
		out[0] = 0.0;
		out[1] = 0.0;
		out[2] = l1(x)*l0(y);
		break;
	case 10:
		out[0] = 0.0;
		out[1] = 0.0;
		out[2] = l0(x)*l1(y);
		break;
	case 11:
		out[0] = 0.0;
		out[1] = 0.0;
		out[2] = l1(x)*l1(y);
	}
}
//---------------------------------------------------------------------------
//-- Базисная функция на векторном параллелепипеде, взятая с весом.
//-- Параллелепипед - любой (не шаблонный), но координаты точки - локальные ([-1,1]^3)
//---------------------------------------------------------------------------
void T_Brick::Basis_func_on_vec_par(long n, double ves, double *in, double *out)
{
	// сначала вычисляем значение базисной функции на шаблонном эл-те
	Basis_func_on_reference_vec_par(n, in, out);

	// домножаем на вес и диагональный эл-т матрицы Якоби
	// (т.к. параллелепипед => всю матрицу Якоби вычислять не нужно)
	switch(n)
	{
	case 0:
	case 1:
	case 2:
	case 3:
		out[0] *= 2.0/hx*ves;
		break;
	case 4:
	case 5:
	case 6:
	case 7:
		out[1] *= 2.0/hy*ves;
		break;
	case 8:
	case 9:
	case 10:
	case 11:
		out[2] *= 2.0/hz*ves;
	}
}
//---------------------------------------------------------------------------
//-- Базисная функция на векторном шестиграннике, взятая с весом.
//-- Шестигранник - любой (не шаблонный), но координаты точки - локальные ([-1,1]^3)
//---------------------------------------------------------------------------
void T_Brick::Basis_func_on_vec_hex(long n, double ves, double *in, double *out)
{
	if(type_of_hex > 30) // векторный шестигранник
	{
		double temp[3];

		Basis_func_on_reference_vec_par(n, in, temp); // на шаблонном эл-те

		Calc_J(in[0], in[1], in[2]);
		Mult_Plot((double*)J_1_T, temp, out, 3); // домножаем на J^{-1}
	}
	else // векторный параллелепипед (матрицу Якоби вычислять не нужно)
	{
		Basis_func_on_vec_par(n, ves, in, out);
	}
}
//---------------------------------------------------------------------------
//---- Вычисляет значение векторного поля внутри произвольного шестигранника
//---- (суммируются все 12 базисных функций с весами).
//---- Координаты точки - локальные
//---------------------------------------------------------------------------
void T_Brick::Calc_value_inside_hex(double *ves, double *in, double *out)
{
	long i, j;
	double temp[3], temp2[3];

	for(i=0; i<3; i++)
		out[i] = 0.0;

	if(type_of_hex > 30) // векторный шестигранник
	{
		Calc_J(in[0], in[1], in[2]); // матрицу Якоби вычисляем один раз

		for (i=0; i<12; i++)
		{
			Basis_func_on_reference_vec_par(i, in, temp);

			Mult_Plot((double*)J_1_T, temp, temp2, 3);
			for (j=0; j<3; j++)
				out[j] += temp2[j]*ves[i];
		}
	} 
	else // векторный параллелепипед
	{
		for (i=0; i<12; i++)
		{
			Basis_func_on_vec_par(i, ves[i], in, temp);
			for (j=0; j<3; j++)
				out[j] += temp[j];
		}
	}
}
//---------------------------------------------------------------------------
//---- Вычисляет значение ротора векторного поля внутри произвольного шестигранника
//---- (суммируются все 12 роторов с весами).
//---- Координаты точки - локальные
//---------------------------------------------------------------------------
void T_Brick::Calc_rotor_inside_hex(double *ves, double *in, double *out)
{
	long i, j, k;
	double temp[3], temp2[3];

	for(i=0; i<3; i++)
		out[i] = 0.0;

	if(type_of_hex > 30) // векторный шестигранник
	{
		Calc_J(in[0], in[1], in[2]); // матрицу Якоби вычисляем один раз

		for (i=0; i<12; i++)
		{
			temp[0] = temp[1] = temp[2] = 0.0;
			Rot_of_basis_func_on_reference_vec_par(i, in, temp2);
			for(k=0; k<3; k++)
				temp[k] +=temp2[k]*ves[i];

			Mult_Plot((double*)J, temp, temp2, 3);
			for (j=0; j<3; j++)
				out[j] += temp2[j]/det_J_abs;
		}
	}
	else // векторный параллелепипед
	{
		for (i=0; i<12; i++)
		{
			Rot_of_basis_func_on_vec_par(i, ves[i], in, temp);

 			switch(i)
 			{
 			case 0:
 			case 1:
 			case 2:
 			case 3:
 				temp[1] *= 2.0/hx;
 				temp[2] *= 2.0/hx;
 				break;
 			case 4:
 			case 5:
 			case 6:
 			case 7:
 				temp[0] *= 2.0/hy;
 				temp[2] *= 2.0/hy;
 				break;
 			case 8:
 			case 9:
 			case 10:
 			case 11:
 				temp[0] *= 2.0/hz;
 				temp[1] *= 2.0/hz;
 			}

			for (j=0; j<3; j++)
				out[j] += temp[j];
		}
	}	
}
//---------------------------------------------------------------------------------------
//---- Вычислить локальную матрицу массы на шестиграннике численно через Гаусс-3 --------
//---------------------------------------------------------------------------------------
void T_Brick::Calc_local_matrix_c_for_hexahedron()
{
	long i, j, k;
	double gauss_3_mult;
	double sig_phi[3],sig_phi0[3];

	if (type_of_hex > 30) // векторный шестигранник
	{
		// обнуляем
		for(i=0; i<12; i++)
		for(j=0; j<=i; j++)
			c[i][j] = 0.0;

		for(k=0; k<27; k++) // по числу точек интегрирования
		{
			Calc_J(k); // вычисляем матрицу Якоби
			gauss_3_mult = gauss_3_A_all[k]*det_J_abs; // A_i*A_j*A_k*|J|

			for(j=0; j<12; j++) // базисные ф-ции преобразуются по правилу...
				Mult_Plot((double*)J_1_T,
				(double*)gauss_3_phi_hex_vec[k][j], (double*)phi_all[j], 3);			

			for(i=0; i<12; i++)
			for(j=0; j<=i; j++)
				c[i][j] += Scal((double*)phi_all[i], (double*)phi_all[j], 3)*gauss_3_mult;
		}

		// матрица симметричная, эл-ты верхнего и нижнего треугольника совпадают
		for(i=0; i<12; i++)
		for(j=i+1; j<12; j++)
			c[i][j] = c[j][i];

		// обнуляем компоненты тензора анизотропии
		for(i=0; i<12; i++)
			for(j=0; j<=i; j++)
			{
				cSigma[i][j] = 0.0;
				cSigma0[i][j] = 0.0;
			}

		for(k=0; k<27; k++) // по числу точек интегрирования
		{
			Calc_J(k); // вычисляем матрицу Якоби
			gauss_3_mult = gauss_3_A_all[k]*det_J_abs; // A_i*A_j*A_k*|J|

			for(j=0; j<12; j++) // базисные ф-ции преобразуются по правилу...
				Mult_Plot((double*)J_1_T,
				(double*)gauss_3_phi_hex_vec[k][j], (double*)phi_all[j], 3);			

			for(i=0; i<12; i++)
			{
				// умножение базисных функций на компоненты тензора анизотропии

				sig_phi[0] = Scal(sigmaTensor.val[0], phi_all[i], 3);
				sig_phi[1] = Scal(sigmaTensor.val[1], phi_all[i], 3);
				sig_phi[2] = Scal(sigmaTensor.val[2], phi_all[i], 3);

				sig_phi0[0] = Scal(s_s0.val[0], phi_all[i], 3);
				sig_phi0[1] = Scal(s_s0.val[1], phi_all[i], 3);
				sig_phi0[2] = Scal(s_s0.val[2], phi_all[i], 3);

				for(j=0; j<=i; j++)
				{
					cSigma[i][j] += Scal(sig_phi, phi_all[j], 3)*gauss_3_mult;
					cSigma0[i][j] += Scal(sig_phi0, phi_all[j], 3)*gauss_3_mult;
				}
			}
		}

		// матрица симметричная, эл-ты верхнего и нижнего треугольника совпадают
		for(i=0; i<12; i++)
			for(j=i+1; j<12; j++)
			{
				cSigma[i][j] = cSigma[j][i];
				cSigma0[i][j] = cSigma0[j][i];
			}
	} 
	else // векторный параллелепипед
	{
		Compute_Local_Matrix_C();
	}
}
//---------------------------------------------------------------------------------------
//---- Вычислить локальную матрицу жёсткости на шестиграннике численно через Гаусс-3 ----
//---------------------------------------------------------------------------------------
void T_Brick::Calc_local_matrix_b_for_hexahedron()
{
	long i, j, k;
	double gauss_3_mult;

	if (type_of_hex > 30) // векторный шестигранник
	{
		// обнуляем
		for(i=0; i<12; i++)
		for(j=0; j<=i; j++)
			b[i][j] = 0.0;

		for(k=0; k<27; k++) // по числу точек интегрирования
		{
			Calc_J(k); // вычисляем матрицу Якоби
			gauss_3_mult = gauss_3_A_all[k]/det_J_abs; // A_i*A_j*A_k/|J|

			for(j=0; j<12; j++) // роторы от базисных ф-ций преобразуются по правилу...
				Mult_Plot((double*)J,
				(double*)gauss_3_rot_phi_hex_vec[k][j], (double*)rot_all[j], 3);			

			for(i=0; i<12; i++)
			for(j=0; j<=i; j++)
				b[i][j] += Scal((double*)rot_all[i], (double*)rot_all[j], 3)*gauss_3_mult;
		}

		// матрица симметричная, эл-ты верхнего и нижнего треугольника совпадают
		for(i=0; i<12; i++) 
		for(j=i+1; j<12; j++)
			b[i][j] = b[j][i];
	} 
	else // векторный параллелепипед
	{
		Compute_Local_Matrix_B();
	}
}
//---------------------------------------------------------------------------------
//--- Вычисляет ротор от базисной функции на параллелепипеде с весом.
//--- Координаты точки - локальные
//---------------------------------------------------------------------------------
void T_Brick::Rot_of_basis_func_on_vec_par(long n, double ves, double *in, double *out)
{
	double x, y, z;
	double in_global[3];
	
	Mapping(in, in_global); // получаем глобальные координаты из локальных

	x = in_global[0];
	y = in_global[1];
	z = in_global[2];

	switch(n)
	{
	case 0:
		out[0] = 0.0;
		out[1] = -(1.0/2.0-(2.0*y-2.0*yk-hy)/hy/2.0)/hz*ves;
		out[2] = 1/hy*(1.0/2.0-(2.0*z-2.0*zk-hz)/hz/2.0)*ves;
		break;
	case 1:
		out[0] = 0.0;
		out[1] = -((2.0*y-2.0*yk-hy)/hy/2.0+1.0/2.0)/hz*ves;
		out[2] = -1/hy*(1.0/2.0-(2.0*z-2.0*zk-hz)/hz/2.0)*ves;
		break;
	case 2:
		out[0] = 0.0;
		out[1] = (1.0/2.0-(2.0*y-2.0*yk-hy)/hy/2.0)/hz*ves;
		out[2] = 1/hy*((2.0*z-2.0*zk-hz)/hz/2.0+1.0/2.0)*ves;
		break;
	case 3:
		out[0] = 0.0;
		out[1] = ((2.0*y-2.0*yk-hy)/hy/2.0+1.0/2.0)/hz*ves;
		out[2] = -1/hy*((2.0*z-2.0*zk-hz)/hz/2.0+1.0/2.0)*ves;
		break;
	case 4:
		out[0] = 1/hz*(1.0/2.0-(2.0*x-2.0*xk-hx)/hx/2.0)*ves;
		out[1] = 0.0;
		out[2] = -(1.0/2.0-(2.0*z-2.0*zk-hz)/hz/2.0)/hx*ves;
		break;
	case 5:
		out[0] = -1/hz*(1.0/2.0-(2.0*x-2.0*xk-hx)/hx/2.0)*ves;
		out[1] = 0.0;
		out[2] = -((2.0*z-2.0*zk-hz)/hz/2.0+1.0/2.0)/hx*ves;
		break;
	case 6:
		out[0] = 1/hz*((2.0*x-2.0*xk-hx)/hx/2.0+1.0/2.0)*ves;
		out[1] = 0.0;
		out[2] = (1.0/2.0-(2.0*z-2.0*zk-hz)/hz/2.0)/hx*ves;
		break;
	case 7:
		out[0] = -1/hz*((2.0*x-2.0*xk-hx)/hx/2.0+1.0/2.0)*ves;
		out[1] = 0.0;
		out[2] = ((2.0*z-2.0*zk-hz)/hz/2.0+1.0/2.0)/hx*ves;
		break;
	case 8:
		out[0] = -(1.0/2.0-(2.0*x-2.0*xk-hx)/hx/2.0)/hy*ves;
		out[1] = 1/hx*(1.0/2.0-(2.0*y-2.0*yk-hy)/hy/2.0)*ves;
		out[2] = 0.0;
		break;
	case 9:
		out[0] = -((2.0*x-2.0*xk-hx)/hx/2.0+1.0/2.0)/hy*ves;
		out[1] = -1/hx*(1.0/2.0-(2.0*y-2.0*yk-hy)/hy/2.0)*ves;
		out[2] = 0.0;
		break;
	case 10:
		out[0] = (1.0/2.0-(2.0*x-2.0*xk-hx)/hx/2.0)/hy*ves;
		out[1] = 1/hx*((2.0*y-2.0*yk-hy)/hy/2.0+1.0/2.0)*ves;
		out[2] = 0.0;
		break;
	case 11:
		out[0] = ((2.0*x-2.0*xk-hx)/hx/2.0+1.0/2.0)/hy*ves;
		out[1] = -1/hx*((2.0*y-2.0*yk-hy)/hy/2.0+1.0/2.0)*ves;
		out[2] = 0.0;
	}
}
//----------------------------------------------------------------------------------------
//- Вычисляет ротор от векторной базисной функции на шаблонном параллелепипеде [-1,1]^3
//- Координаты точки - локальные
//----------------------------------------------------------------------------------------
void T_Brick::Rot_of_basis_func_on_reference_vec_par(long n, double *in, double *out)
{
	double x, y, z;

	x = in[0];
	y = in[1];
	z = in[2];

	switch(n)
	{
	case 0:
		out[0] = 0.0;
		out[1] = -1.0/4.0+y/4.0;
		out[2] = 1.0/4.0-z/4.0;
		break;
	case 1:
		out[0] = 0.0;
		out[1] = -y/4.0-1.0/4.0;
		out[2] = -1.0/4.0+z/4.0;
		break;
	case 2:
		out[0] = 0.0;
		out[1] = 1.0/4.0-y/4.0;
		out[2] = 1.0/4.0+z/4.0;
		break;
	case 3:
		out[0] = 0.0;
		out[1] = y/4.0+1.0/4.0;
		out[2] = -z/4.0-1.0/4.0;
		break;
	case 4:
		out[0] = 1.0/4.0-x/4.0;
		out[1] = 0.0;
		out[2] = -1.0/4.0+z/4.0;
		break;
	case 5:
		out[0] = -1.0/4.0+x/4.0;
		out[1] = 0.0;
		out[2] = -z/4.0-1.0/4.0;
		break;
	case 6:
		out[0] = 1.0/4.0+x/4.0;
		out[1] = 0.0;
		out[2] = 1.0/4.0-z/4.0;
		break;
	case 7:
		out[0] = -x/4.0-1.0/4.0;
		out[1] = 0.0;
		out[2] = 1.0/4.0+z/4.0;
		break;
	case 8:
		out[0] = -1.0/4.0+x/4.0;
		out[1] = 1.0/4.0-y/4.0;
		out[2] = 0.0;
		break;
	case 9:
		out[0] = -x/4.0-1.0/4.0;
		out[1] = -1.0/4.0+y/4.0;
		out[2] = 0.0;
		break;
	case 10:
		out[0] = 1.0/4.0-x/4.0;
		out[1] = y/4.0+1.0/4.0;
		out[2] = 0.0;
		break;
	case 11:
		out[0] = 1.0/4.0+x/4.0;
		out[1] = -y/4.0-1.0/4.0;
		out[2] = 0.0;
	}
}
//------------------------------------------------------------------------
// только x-компонента ротора базисных ф-й на шаблонном параллелепипеде [-1,1]^3
// Координаты точки - локальные
//------------------------------------------------------------------------
void T_Brick::Rotx_of_basis_func_on_reference_vec_par(long n, double x, double *out)
{
	switch(n)
	{
	case 0:
		*out = 0;
		break;
	case 1:
		*out = 0;
		break;
	case 2:
		*out = 0;
		break;
	case 3:
		*out = 0;
		break;
	case 4:
		*out = 1.0/4.0-x/4.0;
		break;
	case 5:
		*out = -1.0/4.0+x/4.0;
		break;
	case 6:
		*out = 1.0/4.0+x/4.0;
		break;
	case 7:
		*out = -x/4.0-1.0/4.0;
		break;
	case 8:
		*out = -1.0/4.0+x/4.0;
		break;
	case 9:
		*out = -x/4.0-1.0/4.0;
		break;
	case 10:
		*out = 1.0/4.0-x/4.0;
		break;
	case 11:
		*out = 1.0/4.0+x/4.0;
	}
}
//------------------------------------------------------------------------
// только y-компонента ротора базисных ф-й на шаблонном параллелепипеде [-1,1]^3
// Координаты точки - локальные
//------------------------------------------------------------------------
void T_Brick::Roty_of_basis_func_on_reference_vec_par(long n, double y, double *out)
{
	switch(n)
	{
	case 0:
		*out = -1.0/4.0+y/4.0;
		break;
	case 1:
		*out = -y/4.0-1.0/4.0;
		break;
	case 2:
		*out = 1.0/4.0-y/4.0;
		break;
	case 3:
		*out = y/4.0+1.0/4.0;
		break;
	case 4:
		*out = 0;
		break;
	case 5:
		*out = 0;
		break;
	case 6:
		*out = 0;
		break;
	case 7:
		*out = 0;
		break;
	case 8:
		*out = 1.0/4.0-y/4.0;
		break;
	case 9:
		*out = -1.0/4.0+y/4.0;
		break;
	case 10:
		*out = y/4.0+1.0/4.0;
		break;
	case 11:
		*out = -y/4.0-1.0/4.0;
	}
}
//------------------------------------------------------------------------
// только z-компонента ротора базисных ф-й на шаблонном параллелепипеде [-1,1]^3
// Координаты точки - локальные
//------------------------------------------------------------------------
void T_Brick::Rotz_of_basis_func_on_reference_vec_par(long n, double z, double *out)
{
	switch(n)
	{
	case 0:
		*out = 0.25-z*0.25;
		break;
	case 1:
		*out = -0.25+z*0.25;
		break;
	case 2:
		*out = 0.25+z*0.25;
		break;
	case 3:
		*out = -z*0.25-0.25;
		break;
	case 4:
		*out = -0.25+z*0.25;
		break;
	case 5:
		*out = -z*0.25-0.25;
		break;
	case 6:
		*out = 0.25-z*0.25;
		break;
	case 7:
		*out = 0.25+z*0.25;
		break;
	case 8:
		*out = 0.0;
		break;
	case 9:
		*out = 0.0;
		break;
	case 10:
		*out = 0.0;
		break;
	case 11:
		*out = 0.0;
	}
}
//-----------------------------------------------------------------------------
//----- преобразование шестигранника из шаблонного в произвольный
//----- (точки шестигранника преобразуются по этому правилу)
//----- т.е. на входе - локальные координаты точки, на выходе - глобальные
//-----------------------------------------------------------------------------
void T_Brick::Mapping(double *in, double *out)
{
	long i;
	double temp;

	out[0] = out[1] = out[2] = 0.0;

	for(i=0; i<8; i++)
	{
		temp = Phi_node(i, in[0], in[1], in[2]);
		out[0] += x[i]*temp;
		out[1] += y[i]*temp;
		out[2] += z[i]*temp;
	}
}
//-----------------------------------------------------------------------------
// вычисляет матрицу Якоби в случае параллелепипеда
//-----------------------------------------------------------------------------
void T_Brick::Calc_J_in_parallelepiped()
{
	// на параллелепипеде эл-ты матрицы Якоби постоянны,
	// поэтому всё равно в какой точке вычислять

	// эл-ты матрицы Якоби
	J[0][0] = hx*0.5;
	J[1][1] = hy*0.5;
	J[2][2] = hz*0.5;

	J[0][1] = J[0][2] = J[1][0] = J[1][2] = J[2][0] = J[2][1] = 0.0;

	// вычисляем Якобиан (определитель)
	det_J = hx*hy*hz/8.0;

	// модуль Якобиана
	det_J_abs = fabs(det_J);

	// матрица, обратная к матрице Якоби (и транспонированная)
	J_1_T[0][0] = J_1[0][0] = 2.0/hx;
	J_1_T[1][1] = J_1[1][1] = 2.0/hy;
	J_1_T[2][2] = J_1[2][2] = 2.0/hz;
	J_1_T[0][1] = J_1_T[0][2] = J_1_T[1][0] = J_1_T[1][2] = J_1_T[2][0] = J_1_T[2][1] =
		J_1[0][1] = J_1[0][2] = J_1[1][0] = J_1[1][2] = J_1[2][0] = J_1[2][1] =	0.0;
}
//-----------------------------------------------------------------------------
// выдать z-компоненту ротора в точках Гаусса в верхней грани шестигранника
// сразу для трёх решений (для 3-слойной схемы по времени)
//-----------------------------------------------------------------------------
void T_Brick::Get_rotz_on_face(double *ves1, double *ves2, double *ves3,
							   double *out1, double *out2, double *out3)
{
	if(type_of_hex > 30) // векторный шестигранник
	{
		double in[3];
		double temp2[3];
		double t;
		int npoint;
		long i;

		in[2] = 1.0;
		for(npoint=0; npoint<9; npoint++)
		{
			in[0] = gauss_3_t_2d[npoint][0];
			in[1] = gauss_3_t_2d[npoint][1];

			Calc_J(in[0], in[1], in[2]); // матрицу Якоби вычисляем один раз для каждой точки Гаусса

			out1[npoint] = out2[npoint] = out3[npoint] = 0.0;

			for (i=0; i<12; i++)
			{
				Rot_of_basis_func_on_reference_vec_par(i, in, temp2);

				// так как требуется только z-компонента ротора,
				// умножается только последняя строчка матрицы Якоби
				t = (J[2][0]*temp2[0] + J[2][1]*temp2[1] + J[2][2]*temp2[2])/det_J_abs;
				out1[npoint] += t*ves1[i];
				out2[npoint] += t*ves2[i];
				out3[npoint] += t*ves3[i];
			}
		}
	}
	else // параллелепипед 
	{
		double t = 2.0/(hx*hy);

		// rot_z на верхней грани=const и определяется 4-мя базисными
		// функциями, ассоциированными с рёбрами, лежащими в этой грани
		out1[0]=out1[1]=out1[2]=out1[3]=out1[4]=out1[5]=
			out1[6]=out1[7]=out1[8] =
				t*ves1[2] - t*ves1[3] - t*ves1[5] + t*ves1[7];

		out2[0]=out2[1]=out2[2]=out2[3]=out2[4]=out2[5]=
			out2[6]=out2[7]=out2[8] =
			t*ves2[2] - t*ves2[3] - t*ves2[5] + t*ves2[7];

		out3[0]=out3[1]=out3[2]=out3[3]=out3[4]=out3[5]=
			out3[6]=out3[7]=out3[8] =
			t*ves3[2] - t*ves3[3] - t*ves3[5] + t*ves3[7];
	}
}
//-----------------------------------------------------------------------------
// выдать компоненты ротора в точках Гаусса в верхней грани шестигранника
// сразу для трёх решений (для 3-слойной схемы по времени)
//-----------------------------------------------------------------------------
void T_Brick::Get_rot_on_face(double *ves1, double *ves2, double *ves3,
							  double *out1, double *out2, double *out3,
							  double *outx1, double *outx2, double *outx3,
							  double *outy1, double *outy2, double *outy3)
{
	if(type_of_hex > 30) // векторный шестигранник
	{
		double in[3];
		double temp2[3];
		double t;
		int npoint;
		long i;

		in[2] = 1.0;
		for(npoint=0; npoint<9; npoint++)
		{
			in[0] = gauss_3_t_2d[npoint][0];
			in[1] = gauss_3_t_2d[npoint][1];

			Calc_J(in[0], in[1], in[2]); // матрицу Якоби вычисляем один раз для каждой точки Гаусса

			out1[npoint] = out2[npoint] = out3[npoint] = 0.0;
			outx1[npoint] = outx2[npoint] = outx3[npoint] = 0.0;
			outy1[npoint] = outy2[npoint] = outy3[npoint] = 0.0;

			for (i=0; i<12; i++)
			{
				Rot_of_basis_func_on_reference_vec_par(i, in, temp2);

				// так как требуется только z-компонента ротора,
				// умножается только последняя строчка матрицы Якоби
				t = (J[2][0]*temp2[0] + J[2][1]*temp2[1] + J[2][2]*temp2[2])/det_J_abs;
				out1[npoint] += t*ves1[i];
				out2[npoint] += t*ves2[i];
				out3[npoint] += t*ves3[i];
				
				// для rot_x
				t = (J[0][0]*temp2[0] + J[0][1]*temp2[1] + J[0][2]*temp2[2])/det_J_abs;
				outx1[npoint] += t*ves1[i];
				outx2[npoint] += t*ves2[i];
				outx3[npoint] += t*ves3[i];
				// для rot_y
				t = (J[1][0]*temp2[0] + J[1][1]*temp2[1] + J[1][2]*temp2[2])/det_J_abs;
				outy1[npoint] += t*ves1[i];
				outy2[npoint] += t*ves2[i];
				outy3[npoint] += t*ves3[i];
			}
		}
	}
	else // параллелепипед 
	{
		int i;

		double t = 2.0/(hx*hy);

		// rot_z на верхней грани=const и определяется 4-мя базисными
		// функциями, ассоциированными с рёбрами, лежащими в этой грани
		out1[0]=out1[1]=out1[2]=out1[3]=out1[4]=out1[5]=
			out1[6]=out1[7]=out1[8] =
				t*ves1[2] - t*ves1[3] - t*ves1[5] + t*ves1[7];

		out2[0]=out2[1]=out2[2]=out2[3]=out2[4]=out2[5]=
			out2[6]=out2[7]=out2[8] =
			t*ves2[2] - t*ves2[3] - t*ves2[5] + t*ves2[7];

		out3[0]=out3[1]=out3[2]=out3[3]=out3[4]=out3[5]=
			out3[6]=out3[7]=out3[8] =
			t*ves3[2] - t*ves3[3] - t*ves3[5] + t*ves3[7];

		
		// для rot_x
		for (i=0; i<9; i++)
		{
			RotXOnPar(gauss_3_t_2d[i][0], ves1, &outx1[i], true);
			RotXOnPar(gauss_3_t_2d[i][0], ves2, &outx2[i], true);
			RotXOnPar(gauss_3_t_2d[i][0], ves3, &outx3[i], true);
		}
		// для rot_y
		for (i=0; i<9; i++)
		{
			RotYOnPar(gauss_3_t_2d[i][1], ves1, &outy1[i], true);
			RotYOnPar(gauss_3_t_2d[i][1], ves2, &outy2[i], true);
			RotYOnPar(gauss_3_t_2d[i][1], ves3, &outy3[i], true);
		}
	}
}
//------------------------------------------------------------------------
// z - компонента родля 3-слойной схемы в произвольной точке
// внутри параллелелепипеда
// z-координата - глобальная
//------------------------------------------------------------------------
void T_Brick::RotZOnPar3(double z,
						double *ves1, double *ves2, double *ves3,
						double *out1, double *out2, double *out3)
{
	int i;
	double t, t2;
	double zz;
	double temp;
	
	*out1 = *out2 = *out3 = 0.0;

	zz = Zeta(z); // Rot_z зависит только от z
	t2 = 4.0/(hx*hy);
	
	// ненулевая компонента ротора только у первых восьми баз. функций
	for (i=0; i<8; i++)
	{
		Rotz_of_basis_func_on_reference_vec_par(i, zz, &temp);

		// так как требуется только z-компонента ротора,
		// умножается только последняя строчка матрицы Якоби
		t = temp*t2;
		*out1 += t*ves1[i];
		*out2 += t*ves2[i];
		*out3 += t*ves3[i];
	}
}
//------------------------------------------------------------------------
//
//------------------------------------------------------------------------
void T_Brick::RotXOnPar(double x, double *ves, double *out, bool loc_c)
{
	int i;
	double t, t2;
	double xx;
	double temp;
	
	*out = 0.0;

	if (!loc_c)
		xx = Xi(x);
	else
		xx = x;
	t2 = 4.0/(hy*hz);

	for (i=4; i<12; i++)
	{
		Rotx_of_basis_func_on_reference_vec_par(i, xx, &temp);

		t = temp*t2;
		*out += t*ves[i];
	}
}
//------------------------------------------------------------------------
//
//------------------------------------------------------------------------
void T_Brick::RotYOnPar(double y, double *ves, double *out, bool loc_c)
{
	int i;
	double t, t2;
	double yy;
	double temp;
	
	*out = 0.0;

	if (!loc_c)
		yy = Eta(y);
	else
		yy = y;
	t2 = 4.0/(hx*hz);

	for (i=0; i<12; i++)
	{
		Roty_of_basis_func_on_reference_vec_par(i, yy, &temp);

		t = temp*t2;
		*out += t*ves[i];
	}
}
//------------------------------------------------------------------------
//
//------------------------------------------------------------------------
void T_Brick::RotZOnPar(double z, double *ves, double *out)
{
	int i;
	double t, t2;
	double zz;
	double temp;
	
	*out = 0.0;

	zz = Zeta(z);
	t2 = 4.0/(hx*hy);

	// ненулевая компонента ротора только у первых восьми баз. функций
	for (i=0; i<8; i++)
	{
		Rotz_of_basis_func_on_reference_vec_par(i, zz, &temp);

		// так как требуется только z-компонента ротора,
		// умножается только последняя строчка матрицы Якоби
		t = temp*t2;
		*out += t*ves[i];
	}
}
//------------------------------------------------------------------------
// выдать x,y,z-компоненту поля в точках Гаусса в верхней грани шестигранника
// сразу для трёх решений (для 3-слойной схемы по времени)
//------------------------------------------------------------------------
void T_Brick::Get_x_y_on_face(double *ves1, double *ves2, double *ves3,
							  double *outx1, double *outx2, double *outx3,
							  double *outy1, double *outy2, double *outy3,
							  double *outz1, double *outz2, double *outz3)
{
	double tempx[3];
	double tempy[3];
	double tempz[3];
	int npoint;

	for(npoint=0; npoint<9; npoint++)
	{
		// вычисляем поле на шаблонном эл-те

		// чтобы выдать x-компоненту в верхней грани нужны рёбра 2 и 3 (с нуля)
		tempx[0] = gauss_3_phi_vec_face[npoint][2][0]*ves1[2] + gauss_3_phi_vec_face[npoint][3][0]*ves1[3];
		tempx[1] = gauss_3_phi_vec_face[npoint][2][0]*ves2[2] + gauss_3_phi_vec_face[npoint][3][0]*ves2[3];
		tempx[2] = gauss_3_phi_vec_face[npoint][2][0]*ves3[2] + gauss_3_phi_vec_face[npoint][3][0]*ves3[3];

		// чтобы выдать y-компоненту в верхней грани нужны рёбра 5 и 7 (с нуля)
		tempy[0] = gauss_3_phi_vec_face[npoint][5][1]*ves1[5] + gauss_3_phi_vec_face[npoint][7][1]*ves1[7];
		tempy[1] = gauss_3_phi_vec_face[npoint][5][1]*ves2[5] + gauss_3_phi_vec_face[npoint][7][1]*ves2[7];
		tempy[2] = gauss_3_phi_vec_face[npoint][5][1]*ves3[5] + gauss_3_phi_vec_face[npoint][7][1]*ves3[7];

		// чтобы выдать z-компоненту в верхней грани нужны рёбра 8,9,10,11 (с нуля)
		tempz[0] = gauss_3_phi_vec_face[npoint][8][2]*ves1[8] + gauss_3_phi_vec_face[npoint][9][2]*ves1[9]
				+ gauss_3_phi_vec_face[npoint][10][2]*ves1[10] + gauss_3_phi_vec_face[npoint][11][2]*ves1[11];

		tempz[1] = gauss_3_phi_vec_face[npoint][8][2]*ves2[8] + gauss_3_phi_vec_face[npoint][9][2]*ves2[9]
				+ gauss_3_phi_vec_face[npoint][10][2]*ves2[10] + gauss_3_phi_vec_face[npoint][11][2]*ves2[11];

		tempz[2] = gauss_3_phi_vec_face[npoint][8][2]*ves3[8] + gauss_3_phi_vec_face[npoint][9][2]*ves3[9]
				+ gauss_3_phi_vec_face[npoint][10][2]*ves3[10] + gauss_3_phi_vec_face[npoint][11][2]*ves3[11];

		// домножаем на Якобиан

		if(type_of_hex > 30) // векторный шестигранник
		{
			// матрицу Якоби вычисляем один раз для каждой точки Гаусса
			Calc_J(gauss_3_t_2d[npoint][0], gauss_3_t_2d[npoint][1], 1.0); 

			// x (умножаем 1-ю строчку матрицы Якоби)
			outx1[npoint] = (J_1_T[0][0]*tempx[0] + J_1_T[0][1]*tempy[0] + J_1_T[0][2]*tempz[0]);
			outx2[npoint] = (J_1_T[0][0]*tempx[1] + J_1_T[0][1]*tempy[1] + J_1_T[0][2]*tempz[1]);
			outx3[npoint] = (J_1_T[0][0]*tempx[2] + J_1_T[0][1]*tempy[2] + J_1_T[0][2]*tempz[2]);

			// y (умножаем 2-ю строчку матрицы Якоби)
			outy1[npoint] = (J_1_T[1][0]*tempx[0] + J_1_T[1][1]*tempy[0] + J_1_T[1][2]*tempz[0]);
			outy2[npoint] = (J_1_T[1][0]*tempx[1] + J_1_T[1][1]*tempy[1] + J_1_T[1][2]*tempz[1]);
			outy3[npoint] = (J_1_T[1][0]*tempx[2] + J_1_T[1][1]*tempy[2] + J_1_T[1][2]*tempz[2]);

			// y (умножаем 3-ю строчку матрицы Якоби)
			outz1[npoint] = (J_1_T[2][0]*tempx[0] + J_1_T[2][1]*tempy[0] + J_1_T[2][2]*tempz[0]);
			outz2[npoint] = (J_1_T[2][0]*tempx[1] + J_1_T[2][1]*tempy[1] + J_1_T[2][2]*tempz[1]);
			outz3[npoint] = (J_1_T[2][0]*tempx[2] + J_1_T[2][1]*tempy[2] + J_1_T[2][2]*tempz[2]);
		}
		else // векторный параллелепипед
		{
			outx1[npoint] = tempx[0]*2.0/hx;
			outx2[npoint] = tempx[1]*2.0/hx;
			outx3[npoint] = tempx[2]*2.0/hx;

			outy1[npoint] = tempy[0]*2.0/hy;
			outy2[npoint] = tempy[1]*2.0/hy;
			outy3[npoint] = tempy[2]*2.0/hy;

			outz1[npoint] = tempz[0]*2.0/hz;
			outz2[npoint] = tempz[1]*2.0/hz;
			outz3[npoint] = tempz[2]*2.0/hz;
		}
	}
}
//--------------------------------------------------------------------------------
//--- Замена переменных для базисных функций, заданных на [-1, 1]
//--- (из глобальных координат --> локальные)
//--------------------------------------------------------------------------------
void T_Brick::Transformation_of_variables(double *in, double *out)
{
	out[0] = Xi(in[0]);
	out[1] = Eta(in[1]);
	out[2] = Zeta(in[2]);
}
//----
void T_Brick::Transformation_of_variables(double *x, double *y, double *z)
{
	x[0] = Xi(x[0]); 
	y[0] = Eta(y[0]); 
	z[0] = Zeta(z[0]); 
}
//----
double T_Brick::Xi(double x)
{
	return (2.0*x - xk - xk1)/hx;
}
//----
double T_Brick::Eta(double y)
{
	return (2.0*y - yk - yk1)/hy;
}
//----
double T_Brick::Zeta(double z)
{
	return (2.0*z - zk - zk1)/hz;
}
//------------------------------------------------------------------------
double T_Brick::GetValueInHexCenter(double *q)
{
	return 0.125*(q[0]+q[1]+q[2]+q[3]+q[4]+q[5]+q[6]+q[7]);
}
//------------------------------------------------------------------------
void T_Brick::GetGradInHexCenter(double *q, double *out, int cj)
{
	int i,j;

	if(cj)Calc_J_Node_2(0.5,0.5,0.5);
	Calc_V_Node_2(q,0.5,0.5,0.5);
	out[0]=out[1]=out[2]=0;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			out[i]+=J_1_T[i][j]*V[j];
		}
	}
}
//------------------------------------------------------------------------
double T_Brick::ScalarFieldOnPar(double x, double y, double z, double *ves)
{
	int i;
	double res = 0.0;
	double coordGlobal[3], coordLocal[3];

	coordGlobal[0] = x;
	coordGlobal[1] = y;
	coordGlobal[2] = z;
	Transformation_of_variables(coordGlobal, coordLocal);

	for (i=0; i<8; i++)
		res += Phi_node(i, coordLocal[0], coordLocal[1], coordLocal[2])*ves[i];

	return res;
}
void T_Brick::ScalarFieldOnParCff(double x, double y, double z, double *cff)
{
	int i;
	double res = 0.0;
	double coordGlobal[3], coordLocal[3];

	coordGlobal[0] = x;
	coordGlobal[1] = y;
	coordGlobal[2] = z;
	Transformation_of_variables(coordGlobal, coordLocal);

	for (i=0; i<8; i++)
		cff[i] = Phi_node(i, coordLocal[0], coordLocal[1], coordLocal[2]);
}
//------------------------------------------------------------------------
double T_Brick::DxOfScalarFieldOnPar(double x, double y, double z, double *ves)
{
	int i;
	double res = 0.0;
	double coordGlobal[3], coordLocal[3];

	coordGlobal[0] = x;
	coordGlobal[1] = y;
	coordGlobal[2] = z;
	Transformation_of_variables(coordGlobal, coordLocal);

	for (i=0; i<8; i++)
		res += dPhi_node(i, 0, coordLocal[0], coordLocal[1], coordLocal[2])*ves[i];
	res = res*2.0/hx;
	
	return res;
}
//------------------------------------------------------------------------
double T_Brick::DyOfScalarFieldOnPar(double x, double y, double z, double *ves)
{
	int i;
	double res = 0.0;
	double coordGlobal[3], coordLocal[3];

	coordGlobal[0] = x;
	coordGlobal[1] = y;
	coordGlobal[2] = z;
	Transformation_of_variables(coordGlobal, coordLocal);

	for (i=0; i<8; i++)
		res += dPhi_node(i, 1, coordLocal[0], coordLocal[1], coordLocal[2])*ves[i];
	res = res*2.0/hy;

	return res;
}
//------------------------------------------------------------------------
double T_Brick::DzOfScalarFieldOnPar(double x, double y, double z, double *ves)
{
	int i;
	double res = 0.0;
	double coordGlobal[3], coordLocal[3];

	coordGlobal[0] = x;
	coordGlobal[1] = y;
	coordGlobal[2] = z;
	Transformation_of_variables(coordGlobal, coordLocal);

	for (i=0; i<8; i++)
		res += dPhi_node(i, 2, coordLocal[0], coordLocal[1], coordLocal[2])*ves[i];
	res = res*2.0/hz;

	return res;
}
//------------------------------------------------------------------------
// выдать значение векторного поля на параллелепипеде
// (координаты точки - глобальные)
//------------------------------------------------------------------------
void T_Brick::VectorFieldOnPar(double x, double y, double z, double *ves,
								 double *x_out, double *y_out, double *z_out)
{
	double out[3], in[3];

	Transformation_of_variables(&x, &y, &z);
	in[0] = x;
	in[1] = y;
	in[2] = z;

	Calc_value_inside_hex(ves, in, out);
	*x_out = out[0];
	*y_out = out[1];
	*z_out = out[2];
}
//------------------------------------------------------------------------
// выдать x-компоненту поля с параллелепипеда сразу для трёх значений весов
//------------------------------------------------------------------------
void T_Brick::VectorFieldXOnPar3(double y, double z, double *ves_j2, double *ves_j1, double *ves_j,
						double *out_j2, double *out_j1, double *out_j)
{
	int i;
	double t[4];

	y = Eta(y);
	z = Zeta(z);

	t[0] = l0(y)*l0(z);
	t[1] = l1(y)*l0(z);
	t[2] = l0(y)*l1(z);
	t[3] = l1(y)*l1(z);

	*out_j2 = 0.0;
	*out_j1 = 0.0;
	*out_j  = 0.0;

	for (i=0; i<4; i++)
	{
		*out_j2 += t[i]*ves_j2[i]; 
		*out_j1 += t[i]*ves_j1[i]; 
		*out_j  += t[i]*ves_j[i]; 
	}

	*out_j2 *= 2.0/hx;
	*out_j1 *= 2.0/hx;
	*out_j  *= 2.0/hx;
}
//------------------------------------------------------------------------
// выдать y-компоненту поля с параллелепипеда сразу для трёх значений весов
//------------------------------------------------------------------------
void T_Brick::VectorFieldYOnPar3(double x, double z, double *ves_j2, double *ves_j1, double *ves_j,
								 double *out_j2, double *out_j1, double *out_j)
{
	int i;
	double t[4];

	x = Xi(x);
	z = Zeta(z);

	t[0] = l0(z)*l0(x);
	t[1] = l1(z)*l0(x);
	t[2] = l0(z)*l1(x);
	t[3] = l1(z)*l1(x);

	out_j2[0] = 0.0;
	out_j1[0] = 0.0;
	out_j[0]  = 0.0;

	for (i=0; i<4; i++)
	{
		out_j2[0] += t[i]*ves_j2[i+4]; 
		out_j1[0] += t[i]*ves_j1[i+4]; 
		out_j[0]  += t[i]*ves_j[i+4]; 
	}

	out_j2[0] *= 2.0/hy;
	out_j1[0] *= 2.0/hy;
	out_j[0]  *= 2.0/hy;
}
//-----------------------------------------------------------------------------
//
//------------------------------------------------------------------------

void T_Brick::Calc_block_local_matrix_and_vector()
{
	long i, j;
	double 
		mult = (sigma0-sigma)*omega,
		mult_dpr = (dpr0-dpr)*omega*omega,
		mult_mu = 1/mu0-1/mu;

	Calc_local_matrix_b_for_hexahedron();
	Calc_local_matrix_c_for_hexahedron();
	Calc_local_vector_for_MT();

	for(i=0; i<12; i++)
	{
		for(j=0; j<12; j++)
		{
			b[i][j] /= mu;
			b[i][j] -= dpr*omega*omega*c0[i][j];
			c[i][j] = cSigma[i][j]*omega;
		}

		if (tasktype!=2)
		{
			g_harm[i*2]  = g_im_sig[i] - g_re[i]*mult_dpr + g_re_b[i]*mult_mu;
			g_harm[i*2+1] = -g_re_sig[i] - g_im[i]*mult_dpr + g_im_b[i]*mult_mu;
		}
		else
		{
			g_harm[i*2]		= g_re_sig[i] - g_im[i]*(dpr-dpr0)*omega;
			g_harm[i*2+1]	= g_im_sig[i] + g_re[i]*(dpr-dpr0)*omega;
		}
	}
}

//-----------------------------------------------------------------------------
//----  вычислить локальный вектор правой части с персчётом 1-мерного поля 
//----   в 3-мерную сетку для МТЗ
//-----------------------------------------------------------------------------
void T_Brick::Calc_local_vector_for_MT()
{
	Calc_asin_acos_at_middle_of_edges();

	Mult_Plot((double*)c0, asin0, g_re, 12);
	Mult_Plot((double*)c0, acos0, g_im, 12);

	Mult_Plot((double*)cSigma0, asin0, g_re_sig, 12);
	Mult_Plot((double*)cSigma0, acos0, g_im_sig, 12);

	Mult_Plot((double*)b, asin0, g_re_b, 12);
	Mult_Plot((double*)b, acos0, g_im_b, 12);
}
//------------------------------------------------------------------------
// вычислить локальный вектор правой части с учётом, что ток=const на элементе
//------------------------------------------------------------------------
void T_Brick::LocalVectHarmConst()
{
	Calc_asin_acos_at_nodes();

	long i, k;
	double gauss_3_mult;
	
	// обнуляем
	for(i=0; i<12; i++)
	{
		g_im[i] = g_re[i] = 0.0;
		g_im_b[i] = g_re_b[i] = 0.0;
	}

	{
		for(k=0; k<27; k++) // по числу точек интегрирования
		{
			Calc_J(k); // вычисляем матрицу Якоби
			gauss_3_mult = gauss_3_A_all[k]*det_J_abs; // A_i*A_j*A_k*|J|

			for(i=0; i<12; i++) // базисные ф-ции преобразуются по правилу...
				Mult_Plot((double*)J_1_T,
				(double*)gauss_3_phi_hex_vec[k][i], (double*)phi_all[i], 3);			

			for(i=0; i<12; i++)
			{
				g_re[i] += Scal(asin0c, (double*)phi_all[i], 3)*gauss_3_mult;
				g_im[i] += Scal(acos0c, (double*)phi_all[i], 3)*gauss_3_mult;
			}
		}
	} 
}
//------------------------------------------------------------------------
// вычислить нормальное поле в вершинах (для гармонических задач)
//------------------------------------------------------------------------
void T_Brick::Calc_asin_acos_at_nodes()
{
	int i, j;
	double coord[3]; // координаты текущей вершины
	double u_sin[3], u_cos[3]; // значение нормального поля в вершине
	double val;

	for (i=0; i<8; i++)
	{
		for (j=0; j<3; j++)
		{
			asin0n[i][j] = 0;
			acos0n[i][j] = 0;
		}
	}

	for (i=0; i<8; i++)
	{
		for (j=0; j<3; j++)
			coord[j] = xyz[nver[num][i]][j];

		{
			if (tasktype==2)
			{
				// значения En есть только в центрах ребер
				u_sin[0]=u_sin[1]=u_sin[2]=u_cos[0]=u_cos[1]=u_cos[2]=0;
			}
			else
			{
			val = Spline(coord[2], n_1d, z_1d, sin_1d);
			u_sin[0] = val*alpha;
			u_sin[1] = val*(1.0 - alpha);
			u_sin[2] = 0.0;

			val = Spline(coord[2], n_1d, z_1d, cos_1d);
			u_cos[0] = val*alpha;
			u_cos[1] = val*(1.0 - alpha);
			u_cos[2] = 0.0;
			}
		}

		for (j=0; j<3; j++)
		{
			asin0n[i][j] += u_sin[j];
			acos0n[i][j] += u_cos[j];
		}
	}

	// суммируем в центр и делим на 8

	for (j=0; j<3; j++)
		asin0c[j] = acos0c[j] = 0.0;

	for (i=0; i<8; i++)
	{
		for (j=0; j<3; j++)
		{
			asin0c[j] += asin0n[i][j];
			acos0c[j] += acos0n[i][j];
		}
	}

	for (j=0; j<3; j++)
	{
		asin0c[j] /= 8;
		acos0c[j] /= 8;
	}

	if (tasktype==2)
	{
		for (j=0; j<3; j++)
			asin0c[j] = acos0c[j] = 0.0;
		// для линии
		// суммируем с ребер в центр и делим на 12
		for (i=0; i<12; i++)
		{
			u_sin[0]=u_sin[1]=u_sin[2]=u_cos[0]=u_cos[1]=u_cos[2]=0;
			for (j=0; j<(int)d->EnForLine[ed[num][i]].size(); j++)
				if (nvkat[num]==d->EnForLine[ed[num][i]][j].mtr)
				{
					const Vec_Prep_Data::EnLine& EnL=d->EnForLine[ed[num][i]+n_edges*ipls][j];
					u_sin[0] = EnL.es[0];
					u_sin[1] = EnL.es[1];
					u_sin[2] = EnL.es[2];
					u_cos[0] = EnL.ec[0];
					u_cos[1] = EnL.ec[1];
					u_cos[2] = EnL.ec[2];
					break;
				}
			for (j=0; j<3; j++)
			{	
				asin0c[j] += u_sin[j];
				acos0c[j] += u_cos[j];
			}
		}
		for (j=0; j<3; j++)
		{
			asin0c[j] /= 12;
			acos0c[j] /= 12;
		}
	}
}

//-----------------------------------------------------------------------------
//-----    вычислить нормальное поле в серединах рёбер (для гармонических задач)
//-----------------------------------------------------------------------------

void normalize(double *v)
{
	double len;
	len=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	v[0]/=len;
	v[1]/=len;
	v[2]/=len;
}

void T_Brick::Calc_asin_acos_at_middle_of_edges()
{
	long i, j;
	double u_sin[3], u_cos[3]; // значение нормального поля в середине ребра
	double x_mid[3]; // координаты точки середины ребра
	double v[3];

	for (i=0; i<12; i++)
	{
		for (j=0; j<3; j++)
		{
			x_mid[j] = 0.5*(xyz[edges[ed[num][i]][0]][j] + xyz[edges[ed[num][i]][1]][j]);
			v[j]=xyz[edges[ed[num][i]][1]][j]-xyz[edges[ed[num][i]][0]][j];
		}
		
		normalize(v);

		{
			if (tasktype==2)
			{
				u_sin[0]=u_sin[1]=u_sin[2]=u_cos[0]=u_cos[1]=u_cos[2]=0;
				for (j=0; j<(int)d->EnForLine[ed[num][i]].size(); j++)
					{
						const Vec_Prep_Data::EnLine& EnL=d->EnForLine[ed[num][i]+n_edges*ipls][j];
						u_sin[0] = EnL.es[0]*v[0];
						u_sin[1] = EnL.es[0]*v[1];
						u_sin[2] = EnL.es[0]*v[2];
						u_cos[0] = EnL.ec[0]*v[0];
						u_cos[1] = EnL.ec[0]*v[1];
						u_cos[2] = EnL.ec[0]*v[2];
						break;
					}
			}
			else
			{
				// 	для реальной задачи МТЗ
  				u_sin[0] = Spline(x_mid[2], n_1d, z_1d, sin_1d)*alpha;
  				u_sin[1] = Spline(x_mid[2], n_1d, z_1d, sin_1d)*(1.0 - alpha);
  				u_sin[2] = 0.0;
		  
  				u_cos[0] = Spline(x_mid[2], n_1d, z_1d, cos_1d)*alpha;
  				u_cos[1] = Spline(x_mid[2], n_1d, z_1d, cos_1d)*(1.0 - alpha);
  				u_cos[2] = 0.0;
			}
		}

		Calc_J(MIDDLE_OF_LOCAL_EDGE[i][0],
			MIDDLE_OF_LOCAL_EDGE[i][1],MIDDLE_OF_LOCAL_EDGE[i][2]);

		this->asin0[i] = Calc_dof((double*)J, u_sin, i);
		this->acos0[i] = Calc_dof((double*)J, u_cos, i);
	}
}
//------------------------------------------------------------------------
void T_Brick::Set_dpr(double dpr)
{
	this->dpr = dpr;
}
//------------------------------------------------------------------------
void T_Brick::Set_dpr0(double dpr0)
{
	this->dpr0 = dpr0;
}
//------------------------------------------------------------------------
void T_Brick::Set_mu0(double mu0)
{
	this->mu0 = mu0;
}
//------------------------------------------------------------------------
void T_Brick::GetVectorFieldNodes(double *ves, double *ax, double *ay, double *az)
{
	int i;
	double a[3];

	for (i=0; i<8; i++)
	{
		Calc_value_inside_hex(ves, (double*)LOCAL_COORDS_OF_NODES[i], a);

		ax[i] = a[0];
		ay[i] = a[1];
		az[i] = a[2];
	}
}
//------------------------------------------------------------------------
