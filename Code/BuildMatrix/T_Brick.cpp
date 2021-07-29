#include "stdafx.h"
#include "T_Brick.h"

extern ofstream logfile;

T_Brick::T_Brick(double *x_coords, double *y_coords, double *z_coords)
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
}

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

	this->alpha = alpha;
	this->n_1d =n_1d;
	this->z_1d = z_1d;
	this->sin_1d = sin_1d;
	this->cos_1d = cos_1d;
}

T_Brick::~T_Brick()
{
}

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
void T_Brick::Compute_Local_Matrix_C()
{
	double mult_x, mult_y, mult_z;
	double t1,t2,t3,t4,t5,t6,t7,t8,t9;

	mult_x = hy*hz/(2.0*hx);
	mult_y = hx*hz/(2.0*hy);
	mult_z = hx*hy/(2.0*hz);

	t1 = 8.0/9.0*mult_x;
	t2 = 4.0/9.0*mult_x;
	t3 = 2.0/9.0*mult_x;
	t4 = 8.0/9.0*mult_y;
	t5 = 4.0/9.0*mult_y;
	t6 = 2.0/9.0*mult_y;
	t7 = 8.0/9.0*mult_z;
	t8 = 4.0/9.0*mult_z;
	t9 = 2.0/9.0*mult_z;

	c[0][0] = t1;
	c[0][1] = t2;
	c[0][2] = t2;
	c[0][3] = t3;
	c[0][4] = 0.0;
	c[0][5] = 0.0;
	c[0][6] = 0.0;
	c[0][7] = 0.0;
	c[0][8] = 0.0;
	c[0][9] = 0.0;
	c[0][10] = 0.0;
	c[0][11] = 0.0;
	c[1][0] = t2;
	c[1][1] = t1;
	c[1][2] = t3;
	c[1][3] = t2;
	c[1][4] = 0.0;
	c[1][5] = 0.0;
	c[1][6] = 0.0;
	c[1][7] = 0.0;
	c[1][8] = 0.0;
	c[1][9] = 0.0;
	c[1][10] = 0.0;
	c[1][11] = 0.0;
	c[2][0] = t2;
	c[2][1] = t3;
	c[2][2] = t1;
	c[2][3] = t2;
	c[2][4] = 0.0;
	c[2][5] = 0.0;
	c[2][6] = 0.0;
	c[2][7] = 0.0;
	c[2][8] = 0.0;
	c[2][9] = 0.0;
	c[2][10] = 0.0;
	c[2][11] = 0.0;
	c[3][0] = t3;
	c[3][1] = t2;
	c[3][2] = t2;
	c[3][3] = t1;
	c[3][4] = 0.0;
	c[3][5] = 0.0;
	c[3][6] = 0.0;
	c[3][7] = 0.0;
	c[3][8] = 0.0;
	c[3][9] = 0.0;
	c[3][10] = 0.0;
	c[3][11] = 0.0;
	c[4][0] = 0.0;
	c[4][1] = 0.0;
	c[4][2] = 0.0;
	c[4][3] = 0.0;
	c[4][4] = t4;
	c[4][5] = t5;
	c[4][6] = t5;
	c[4][7] = t6;
	c[4][8] = 0.0;
	c[4][9] = 0.0;
	c[4][10] = 0.0;
	c[4][11] = 0.0;
	c[5][0] = 0.0;
	c[5][1] = 0.0;
	c[5][2] = 0.0;
	c[5][3] = 0.0;
	c[5][4] = t5;
	c[5][5] = t4;
	c[5][6] = t6;
	c[5][7] = t5;
	c[5][8] = 0.0;
	c[5][9] = 0.0;
	c[5][10] = 0.0;
	c[5][11] = 0.0;
	c[6][0] = 0.0;
	c[6][1] = 0.0;
	c[6][2] = 0.0;
	c[6][3] = 0.0;
	c[6][4] = t5;
	c[6][5] = t6;
	c[6][6] = t4;
	c[6][7] = t5;
	c[6][8] = 0.0;
	c[6][9] = 0.0;
	c[6][10] = 0.0;
	c[6][11] = 0.0;
	c[7][0] = 0.0;
	c[7][1] = 0.0;
	c[7][2] = 0.0;
	c[7][3] = 0.0;
	c[7][4] = t6;
	c[7][5] = t5;
	c[7][6] = t5;
	c[7][7] = t4;
	c[7][8] = 0.0;
	c[7][9] = 0.0;
	c[7][10] = 0.0;
	c[7][11] = 0.0;
	c[8][0] = 0.0;
	c[8][1] = 0.0;
	c[8][2] = 0.0;
	c[8][3] = 0.0;
	c[8][4] = 0.0;
	c[8][5] = 0.0;
	c[8][6] = 0.0;
	c[8][7] = 0.0;
	c[8][8] = t7;
	c[8][9] = t8;
	c[8][10] = t8;
	c[8][11] = t9;
	c[9][0] = 0.0;
	c[9][1] = 0.0;
	c[9][2] = 0.0;
	c[9][3] = 0.0;
	c[9][4] = 0.0;
	c[9][5] = 0.0;
	c[9][6] = 0.0;
	c[9][7] = 0.0;
	c[9][8] = t8;
	c[9][9] = t7;
	c[9][10] = t9;
	c[9][11] = t8;
	c[10][0] = 0.0;
	c[10][1] = 0.0;
	c[10][2] = 0.0;
	c[10][3] = 0.0;
	c[10][4] = 0.0;
	c[10][5] = 0.0;
	c[10][6] = 0.0;
	c[10][7] = 0.0;
	c[10][8] = t8;
	c[10][9] = t9;
	c[10][10] = t7;
	c[10][11] = t8;
	c[11][0] = 0.0;
	c[11][1] = 0.0;
	c[11][2] = 0.0;
	c[11][3] = 0.0;
	c[11][4] = 0.0;
	c[11][5] = 0.0;
	c[11][6] = 0.0;
	c[11][7] = 0.0;
	c[11][8] = t9;
	c[11][9] = t8;
	c[11][10] = t8;
	c[11][11] = t7;
}

void T_Brick::Compute_Local_Matrix_And_Vector(const long what_compute)
{
	long i,j;

	switch(what_compute)
	{

	case 1: // вычислить матрицу массы и вектор
		Compute_Local_Matrix_C();

		for(i=0; i<12; i++)
		{
			for(j=0; j<12; j++)
				a[i][j] = c[i][j]*sigma;
			g[i] = 0.0;
		}
		break;

	case 2: // если требуетс€ только прибавить матрицу жЄсткости
		Compute_Local_Matrix_B();

		for(i=0; i<12; i++)
		{
			for(j=0; j<12; j++)
				a[i][j] = b[i][j]/mu;
			g[i] = 0.0;
		}
		break;

	case 3: // собрать матрицу массы с коэффициентом диэлектрической проницаемости
		Compute_Local_Matrix_C();

		for(i=0; i<12; i++)
		{
			for(j=0; j<12; j++)
				a[i][j] = c[i][j]*dpr;
			g[i] = 0.0;
		}

		break;
	}
}

void T_Brick::Calc_J(int n_of_point)
{
	Calc_J_in_parallelepiped();		
}

void T_Brick::Calc_J(double xi, double eta, double zeta)
{
	Calc_J_in_parallelepiped();
}

void T_Brick::Calc_J_on_face(int n_of_point)
{
	Calc_J_in_parallelepiped();
}

double T_Brick::l0(double x)
{
	return (1.0 - x)*0.5;
}

double T_Brick::l1(double x)
{
	return (x + 1.0)*0.5;
}

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

void T_Brick::Basis_func_on_vec_par(long n, double ves, double *in, double *out)
{
	Basis_func_on_reference_vec_par(n, in, out);

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

void T_Brick::Basis_func_on_vec_hex(long n, double ves, double *in, double *out)
{
	Basis_func_on_vec_par(n, ves, in, out);
}

void T_Brick::Calc_value_inside_hex(double *ves, double *in, double *out)
{
	long i, j;
	double temp[3], temp2[3];

	for(i=0; i<3; i++)
		out[i] = 0.0;

	for (i=0; i<12; i++)
	{
		Basis_func_on_vec_par(i, ves[i], in, temp);
		for (j=0; j<3; j++)
			out[j] += temp[j];
	}
}

void T_Brick::Calc_rotor_inside_hex(double *ves, double *in, double *out)
{
	long i, j, k;
	double temp[3], temp2[3];

	for(i=0; i<3; i++)
		out[i] = 0.0;

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

void T_Brick::Rot_of_basis_func_on_vec_par(long n, double ves, double *in, double *out)
{
	double x, y, z;
	double in_global[3];
	
	x=0.5*(in[0]+1.0)*hx+this->x[0];
	y=0.5*(in[1]+1.0)*hy+this->y[0];
	z=0.5*(in[2]+1.0)*hz+this->z[0];

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

void T_Brick::Calc_J_in_parallelepiped()
{
	// на параллелепипеде эл-ты матрицы якоби посто€нны,
	// поэтому всЄ равно в какой точке вычисл€ть

	// эл-ты матрицы якоби
	J[0][0] = hx*0.5;
	J[1][1] = hy*0.5;
	J[2][2] = hz*0.5;

	J[0][1] = J[0][2] = J[1][0] = J[1][2] = J[2][0] = J[2][1] = 0.0;

	// вычисл€ем якобиан (определитель)
	det_J = hx*hy*hz/8.0;

	// модуль якобиана
	det_J_abs = fabs(det_J);

	// матрица, обратна€ к матрице якоби (и транспонированна€)
	J_1_T[0][0] = J_1[0][0] = 2.0/hx;
	J_1_T[1][1] = J_1[1][1] = 2.0/hy;
	J_1_T[2][2] = J_1[2][2] = 2.0/hz;
	J_1_T[0][1] = J_1_T[0][2] = J_1_T[1][0] = J_1_T[1][2] = J_1_T[2][0] = J_1_T[2][1] =
		J_1[0][1] = J_1[0][2] = J_1[1][0] = J_1[1][2] = J_1[2][0] = J_1[2][1] =	0.0;
}

void T_Brick::Get_rotz_on_face(double *ves1, double *ves2, double *ves3,
							   double *out1, double *out2, double *out3)
{
	double t = 2.0/(hx*hy);

	// rot_z на верхней грани=const и определ€етс€ 4-м€ базисными
	// функци€ми, ассоциированными с рЄбрами, лежащими в этой грани
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
	
	// ненулева€ компонента ротора только у первых восьми баз. функций
	for (i=0; i<8; i++)
	{
		Rotz_of_basis_func_on_reference_vec_par(i, zz, &temp);

		// так как требуетс€ только z-компонента ротора,
		// умножаетс€ только последн€€ строчка матрицы якоби
		t = temp*t2;
		*out1 += t*ves1[i];
		*out2 += t*ves2[i];
		*out3 += t*ves3[i];
	}
}

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

void T_Brick::RotZOnPar(double z, double *ves, double *out)
{
	int i;
	double t, t2;
	double zz;
	double temp;
	
	*out = 0.0;

	zz = Zeta(z);
	t2 = 4.0/(hx*hy);

	// ненулева€ компонента ротора только у первых восьми баз. функций
	for (i=0; i<8; i++)
	{
		Rotz_of_basis_func_on_reference_vec_par(i, zz, &temp);

		// так как требуетс€ только z-компонента ротора,
		// умножаетс€ только последн€€ строчка матрицы якоби
		t = temp*t2;
		*out += t*ves[i];
	}
}

void T_Brick::Transformation_of_variables(double *in, double *out)
{
	out[0] = Xi(in[0]);
	out[1] = Eta(in[1]);
	out[2] = Zeta(in[2]);
}

void T_Brick::Transformation_of_variables(double *x, double *y, double *z)
{
	x[0] = Xi(x[0]); 
	y[0] = Eta(y[0]); 
	z[0] = Zeta(z[0]); 
}

double T_Brick::Xi(double x)
{
	return (2.0*x - xk - xk1)/hx;
}

double T_Brick::Eta(double y)
{
	return (2.0*y - yk - yk1)/hy;
}

double T_Brick::Zeta(double z)
{
	return (2.0*z - zk - zk1)/hz;
}

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

void T_Brick::Calc_block_local_matrix_and_vector()
{
	long i, j;
	double 
		mult = (sigma0-sigma)*omega,
		mult_dpr = (dpr0-dpr)*omega*omega,
		mult_mu = 1/mu0-1/mu;

	Compute_Local_Matrix_B();
	Compute_Local_Matrix_C();
	Calc_local_vector_for_MT();

	for(i=0; i<12; i++)
	{
		for(j=0; j<12; j++)
		{
			b[i][j] /= mu;
			b[i][j] -= dpr*omega*omega*c[i][j];
			c[i][j] *= sigma*omega;
		}
		g_harm[i*2]  = g_im[i]*(-mult)- g_re[i]*mult_dpr + g_re_b[i]*mult_mu;
		g_harm[i*2+1] = g_re[i]*mult  - g_im[i]*mult_dpr + g_im_b[i]*mult_mu;
	}
}

void T_Brick::Calc_local_vector_for_MT()
{
	Calc_asin_acos_at_middle_of_edges();

	Mult_Plot((double*)c, asin0, g_re, 12);
	Mult_Plot((double*)c, acos0, g_im, 12);

	Mult_Plot((double*)b, asin0, g_re_b, 12);
	Mult_Plot((double*)b, acos0, g_im_b, 12);
}

void T_Brick::Calc_asin_acos_at_nodes()
{
	int i, j;
	double coord[3];
	double u_sin[3], u_cos[3];
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

		val = Spline(coord[2], n_1d, z_1d, sin_1d);
		u_sin[0] = val*alpha;
		u_sin[1] = val*(1.0 - alpha);
		u_sin[2] = 0.0;

		val = Spline(coord[2], n_1d, z_1d, cos_1d);
		u_cos[0] = val*alpha;
		u_cos[1] = val*(1.0 - alpha);
		u_cos[2] = 0.0;

		for (j=0; j<3; j++)
		{
			asin0n[i][j] += u_sin[j];
			acos0n[i][j] += u_cos[j];
		}
	}

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
}

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
	double u_sin[3], u_cos[3]; // значение нормального пол€ в середине ребра
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

		u_sin[0] = Spline(x_mid[2], n_1d, z_1d, sin_1d)*alpha;
		u_sin[1] = Spline(x_mid[2], n_1d, z_1d, sin_1d)*(1.0 - alpha);
		u_sin[2] = 0.0;
  
		u_cos[0] = Spline(x_mid[2], n_1d, z_1d, cos_1d)*alpha;
		u_cos[1] = Spline(x_mid[2], n_1d, z_1d, cos_1d)*(1.0 - alpha);
		u_cos[2] = 0.0;

		Calc_J(MIDDLE_OF_LOCAL_EDGE[i][0],
			MIDDLE_OF_LOCAL_EDGE[i][1],MIDDLE_OF_LOCAL_EDGE[i][2]);

		this->asin0[i] = Calc_dof((double*)J, u_sin, i);
		this->acos0[i] = Calc_dof((double*)J, u_cos, i);
	}	
}

void T_Brick::Set_dpr(double dpr)
{
	this->dpr = dpr;
}

void T_Brick::Set_dpr0(double dpr0)
{
	this->dpr0 = dpr0;
}

void T_Brick::Set_mu0(double mu0)
{
	this->mu0 = mu0;
}

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
