#pragma once
#include "Vec_Prep_Data.h"

class T_Brick
{
public:
	double hx, hy, hz; // размеры эл-та
	double xk, xk1, yk, yk1, zk, zk1; // границы параллелепипеда
    long num; // номер конечного элемента в сетке

	double b[12][12]; // локальная матрица жёсткости
	double c[12][12]; // локальная матрица массы

	double f_re[12];  // значения правой части в серединах рёбер

	double a[12][12]; // локальная матрица эл-та  
	double g[12];     // локальный вектор эл-та
	double g_harm[24];
	double g_re[12], g_im[12];
	double g_re_b[12], g_im_b[12];

	// коэффициенты уравнения 
	double mu;
	double mu0;
	double sigma;
	double sigma0;
	double dpr;
	double dpr0;
	long n_mat; // номер материала элемента

	long (*nver)[14]; // номера узлов конечных элементов (13-узловые шестигранники)
	long (*ed)[25];   // эл-ты перечисленные своими рёбрами + терминальные рёбра + тип эл-та
	long *nvkat;      // номера материалов конечных эл-тов
	long (*edges)[2]; // рёбра, заданные 2-мя вершинами
	double (*xyz)[3]; // координаты узлов

	T_Brick(double *x_coords, double *y_coords, double *z_coords);

	T_Brick(long num, long (*nver)[14], long (*ed)[25], long (*edges)[2], double (*xyz)[3], long *nvkat,
		double *sigma3d, double *sigma0, double *mu3d, double omega,
		long alpha, long n_1d, double *z_1d, double *sin_1d, double *cos_1d);

	T_Brick(long num, long (*nver)[14], double (*xyz)[3]);

	~T_Brick();

	void Compute_Local_Matrix_And_Vector(const long what_compute); 
	void Compute_Local_Matrix_B(); 
	void Compute_Local_Matrix_C(); 

	double omega; // циклическая частота для гармонической задачи
	double asin0[12], acos0[12]; // нормальное поле в серединах рёбер 

	long alpha;
	long n_1d;
	double *z_1d;
	double *sin_1d;
	double *cos_1d;

	void Calc_local_vector_for_MT();

	void Calc_asin_acos_at_middle_of_edges();

	void Calc_block_local_matrix_and_vector();

	double x[8], y[8], z[8]; // координаты вершин шестигранника	
	double J[3][3];          // матрица Якоби
	double J_1[3][3];		 // J^{-1}
	double J_1_T[3][3];		 // J^{-T}
	double det_J;            // Якобиан (определитель матрицы Якоби)
	double det_J_abs;        // модуль Якобиана
	double phi_all[12][3];    // базисные функций в точке интегрирования
	double rot_all[12][3];    // роторы от баз. функций в точке интегрирования

	void Mapping(double *in, double *out);

	void Calc_J(int n_of_point); // вычисляет матрицу Якоби в точке Гаусса 
	void Calc_J(double x, double y, double z); // вычисляет матрицу Якоби в произвольной точке 
	void Calc_J_on_face(int n_of_point); // вычисляет Якобиан в точках, расположенных в грани
	void Calc_J_in_parallelepiped(); // вычисляет матрицу Якоби в случае параллелепипеда

	// Вычисляет значение векторного поля внутри произвольного шестигранника
	void Calc_value_inside_hex(double *ves, double *in, double *out); 

	// значение i-й базисные функции на параллелепипеде 
	void Basis_func_on_vec_par(long i, double *in, double *out);

	// значение i-й базисные функции на параллелепипеде с тангенциальными составляющими 2/hx, 2/hy, 2/hz
	void Basis_func_on_vec_par(long i, double ves, double *in, double *out);

	// значение i-й базисной функции на шаблонном параллелепипеде 
	void Basis_func_on_reference_vec_par(long i, double *in, double *out);	

	// значение i-й базисной функции на шестиграннике 
	void Basis_func_on_vec_hex(long i, double ves, double *in, double *out);

	// выдать ротор в произвольной точке внутри шестигранника
	void Calc_rotor_inside_hex(double *ves, double *in, double *out); 

	// значение ротора i-й базисной функции на шаблонном параллелепипеде 
	void Rot_of_basis_func_on_reference_vec_par(long i, double *in, double *out);

	// значение только x-компоненты ротора i-й базисной функции внутри параллелепипеда
	void Rotx_of_basis_func_on_reference_vec_par(long i, double x, double *out);

	// значение только y-компоненты ротора i-й базисной функции внутри параллелепипеда
	void Roty_of_basis_func_on_reference_vec_par(long i, double y, double *out);

	// значение только z-компоненты ротора i-й базисной функции внутри параллелепипеда
	void Rotz_of_basis_func_on_reference_vec_par(long i, double z, double *out);

	// значение ротора i-й базисной функции внутри параллелепипеда
	void Rot_of_basis_func_on_vec_par(long i, double ves, double *in, double *out);

	// выдать z-компоненту ротора в точках Гаусса в верхней грани шестигранника
	// сразу для трёх решений (для 3-слойной схемы по времени)
	void Get_rotz_on_face(double *ves1, double *ves2, double *ves3,
		double *out1, double *out2, double *out3);

	// Шаблонные одномерные базисные функции на [-1, 1]^3
	double l0(double x);
	double l1(double x);

	// Шаблонные узловые базисные функции на [-1, 1]^3
	double Phi_node(long i, double x, double y, double z);

	// производные от узловых шаблонных базисных ф-ций
	double dPhi_node(long i, long j, double x, double y, double z);

	// выдать значение векторного поля на параллелепипеде
	// (координаты точки - глобальные)
	void VectorFieldOnPar(double x, double y, double z, double *ves,
		double *x_out, double *y_out, double *z_out);
	double ScalarFieldOnPar(double x, double y, double z, double *ves);
	double DxOfScalarFieldOnPar(double x, double y, double z, double *ves);
	double DyOfScalarFieldOnPar(double x, double y, double z, double *ves);
	double DzOfScalarFieldOnPar(double x, double y, double z, double *ves);

	// выдать x-компоненту поля с параллелепипеда сразу для трёх значений весов
	void VectorFieldXOnPar3(double y, double z, double *ves_j2, double *ves_j1, double *ves_j,
		double *out_j2, double *out_j1, double *out_j);

	// выдать y-компоненту поля с параллелепипеда сразу для трёх значений весов
	void VectorFieldYOnPar3(double x, double z, double *ves_j2, double *ves_j1, double *ves_j,
		double *out_j2, double *out_j1, double *out_j);

	void RotXOnPar(double x, double *ves, double *out, bool loc_c=false);
	void RotYOnPar(double y, double *ves, double *out, bool loc_c=false);
	void RotZOnPar(double z, double *ves, double *out);
	void RotZOnPar3(double z,
		double *ves1, double *ves2, double *ves3,
		double *out1, double *out2, double *out3);

	// замена переменных
	void Transformation_of_variables(double *in, double *out);
	void Transformation_of_variables(double *x, double *y, double *z);
	double Xi(double x);
	double Eta(double y);
	double Zeta(double z);

	void GetAonFace(double *ves1, double *ves2, double *ves3,
		double *out1, double *out2, double *out3);

	void Set_dpr(double dpr);
	void Set_dpr0(double dpr0);
	void Set_mu0(double mu0);
	Vec_Prep_Data *d;

	double asin0n[8][3], acos0n[8][3]; // нормальное поле в узлах
	double asin0c[3], acos0c[3]; // нормальное поле в центре элемента

	void Calc_asin_acos_at_nodes();

	void GetVectorFieldNodes(double *ves, double *ax, double *ay, double *az);

	int n_edges;
};
//------------------------------------------------------------------------------
const double MIDDLE_OF_LOCAL_EDGE[12][3] = {
	 0.0, -1.0, -1.0,
	 0.0,  1.0, -1.0,
	 0.0, -1.0,  1.0,
	 0.0,  1.0,  1.0,

	-1.0,  0.0, -1.0,
	-1.0,  0.0,  1.0,
	 1.0,  0.0, -1.0,
	 1.0,  0.0,  1.0,

	-1.0, -1.0,  0.0,
	 1.0, -1.0,  0.0,
	-1.0,  1.0,  0.0,
	 1.0,  1.0,  0.0
};
//-----------------------------------------------------------------------
const double LOCAL_COORDS_OF_NODES[8][3] = 
{
	-1.0, -1.0, -1.0,
	 1.0, -1.0, -1.0,
	-1.0,  1.0, -1.0,
	 1.0,  1.0, -1.0,

	-1.0, -1.0,  1.0,
	 1.0, -1.0,  1.0,
	-1.0,  1.0,  1.0,
	 1.0,  1.0,  1.0
};
//-----------------------------------------------------------------------
const double TANGENT_VECTORS_ON_REFERENCE_CUBE[12][3] = {
	1.0, 0.0, 0.0,
	1.0, 0.0, 0.0,
	1.0, 0.0, 0.0,
	1.0, 0.0, 0.0,

	0.0, 1.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 1.0, 0.0,

	0.0, 0.0, 1.0,
	0.0, 0.0, 1.0,
	0.0, 0.0, 1.0,
	0.0, 0.0, 1.0
};
//-----------------------------------------------------------------------
const long REG_EDGES[12][2]={ 	// какие нетерминальные рёбра есть на элементе
	0,1, 2,3, 4,5, 6,7, 0,2, 4,6, 1,3, 5,7, 0,4, 1,5, 2,6, 3,7 };
//-----------------------------------------------------------------------