#pragma once
#include "ListForLine.h"
#include "Vec_Prep_Data.h"

//------------------------------------------------------------------------------------
/*	Класс T_Brick содержит функции для работы с векторными (edge-elements) шестигранниками
	для нестационарной параболической задачи (петля):
	- вычисление локальных матриц жёсткости и массы
	- вычисление вектора правой части через E_нормальное
	- вычисление матрицы Якоби в точках Гаусса и в произвольной точке внутри шестигранника
	- выдача решения и ротора внутри шестигранника или параллелепипеда

	Узловые базисные функции здесь нужны для вычисления преобразования.
	Шаблонный элемент [-1, 1]^3. (как для узловых, так и для векторных).
	Векторные базисные функции с тангенциальными составляющими 2/hx, 2/hy, 2/hz вдоль рёбер.
*/
//------------------------------------------------------------------------------------
class T_Brick
{
public:
	double hx, hy, hz; // размеры эл-та
	double xk, xk1, yk, yk1, zk, zk1; // границы параллелепипеда
    long num; // номер конечного элемента в сетке

	double b[12][12]; // локальная матрица жёсткости
	double c[12][12]; // локальная матрица массы
	double c0[12][12]; // локальная матрица массы
	double cSigma[12][12]; // локальная матрица массы
	double cSigma0[12][12]; // локальная матрица массы

	double f_re[12];  // значения правой части в серединах рёбер

	double a[12][12]; // локальная матрица эл-та  
	double g[12];     // локальный вектор эл-та
	double g_harm[24];
	double g_re[12], g_im[12];
	double g_re_sig[12], g_im_sig[12];
	double g_re_b[12], g_im_b[12];
	double g8[8]; // для узловых

	// коэффициенты уравнения 
	double mu;
	double mu0;
	double sigma;
    Tensor sigmaTensor, sigmaTensor0, s_s0;
	double sigma0;
	double dpr;
	double dpr0;
	long n_mat; // номер материала элемента

	long (*nver)[14]; // номера узлов конечных элементов (13-узловые шестигранники)
	long (*ed)[25];   // эл-ты перечисленные своими рёбрами + терминальные рёбра + тип эл-та
	long *nvkat;      // номера материалов конечных эл-тов
	long (*edges)[2]; // рёбра, заданные 2-мя вершинами
	double (*xyz)[3]; // координаты узлов

	double *En; // нормальное поле (на рёбрах)
	double (*En_nodes)[3]; // в узлах 

	//------------------start--AV постановка--start-------------------------------------------
	long (*nvetr)[20]; // номера узлов(1-8) и ребер(9-20) в массиве nded
	double MtrMass[20][20];
	double MtrMassRight[20][12];
	double MtrGest[20][20];
	double VctRight[40];
	double g_re_sig_av[20], g_im_sig_av[20];

	// конструктор для МТЗ
	T_Brick(long num, long (*nvetr)[20], long (*nver)[14], long (*ed)[25], long (*edges)[2], double (*xyz)[3], long *nvkat,
		Tensor *sigma3d, double *sigma0, double *mu3d, double omega,
		long alpha, long n_1d, double *z_1d, double *sin_1d, double *cos_1d);
	void Calc_block_local_matrix_and_vector_AV();
	// вычислить локальную матрицу жёсткости на параллелепипеде
	void Compute_Local_Matrix_B_AV(); 
	// вычислить локальную матрицу массы на параллелепипеде
	void Compute_Local_Matrix_C_AV(); 
	// вычислить локальную матрицу 20x12 для формирования правой части
	void Compute_Local_Matrix_Pr_AV(); 
	// вычислить локальный вектор правой части с персчётом 1-мерного поля 
	// в 3-мерную сетку для МТЗ
	void Calc_local_vector_for_MT_AV();
	//------------------end--AV постановка--end-------------------------------------------

	int tasktype;

	// этот конструктор обычно используется для выдачи с элемента
	T_Brick(double *x_coords, double *y_coords, double *z_coords, long type_of_hex);

	// конструктор для нестационарной задачи (для вычисления локальных матриц)
	T_Brick(long num, long (*nver)[14], long (*ed)[25], long (*edges)[2], double (*xyz)[3], long *nvkat,
		double *sigma3d, double *sigma0, double *mu3d, double *En);

	// конструктор для МТЗ
	T_Brick(long num, long (*nver)[14], long (*ed)[25], long (*edges)[2], double (*xyz)[3], long *nvkat,
		double *sigma3d, double *sigma0, double *mu3d, double omega,
		long alpha, long n_1d, double *z_1d, double *sin_1d, double *cos_1d);

	T_Brick(long num, long (*nver)[14], double (*xyz)[3]);

	~T_Brick();

	//	вычислить локальную матрицу для нестационарной задачи
	// (матрица массы и вектор или только прибавить матрицу жёсткости)
	void Compute_Local_Matrix_And_Vector(const long what_compute); 

	// вычислить локальную матрицу жёсткости на параллелепипеде
	void Compute_Local_Matrix_B(); 

	// вычислить локальную матрицу массы на параллелепипеде
	void Compute_Local_Matrix_C(); 

	// вычисление вектора правой части через E_нормальное
	void Compute_Local_Vector_For_Anomal_Problem();

	// вычисление вектора правой части через нормальное поле, в случае, когда все коэф-ты разрывны
	void ComputeLocalVectorMuEpsSigma(double *An, double *d2An);

	//	для МТЗ	 

	double omega; // циклическая частота для гармонической задачи
	double asin0[12], acos0[12]; // нормальное поле в серединах рёбер 

	// для одномерки в МТЗ 
	long alpha;
	long n_1d;
	double *z_1d;
	double *sin_1d;
	double *cos_1d;

	// вычислить локальный вектор правой части с персчётом 1-мерного поля 
	// в 3-мерную сетку для МТЗ
	void Calc_local_vector_for_MT();

	// вычислить нормальное поле в серединах рёбер (для МТЗ)
	void Calc_asin_acos_at_middle_of_edges();

	void Calc_block_local_matrix_and_vector();

	void Calc_J_Node_2(double x, double y, double z);
	double dPhi_node_2(long i, long j, double x, double y, double z);
	void Calc_V_Node_2(double *q,double x, double y, double z);
	double V[3];

	// для векторных шестигранников

	double x[8], y[8], z[8]; // координаты вершин шестигранника	
	double J[3][3];          // матрица Якоби
	double J_1[3][3];		 // J^{-1}
	double J_1_T[3][3];		 // J^{-T}
	double det_J;            // Якобиан (определитель матрицы Якоби)
	double det_J_abs;        // модуль Якобиана
	long type_of_hex;        // 0..30 - параллелепипед, 31..61 - шестигранник
	double phi_all[12][3];    // базисные функций в точке интегрирования
	double rot_all[12][3];    // роторы от баз. функций в точке интегрирования

	// преобразование шестигранника из шаблонного в произвольный
	// (получаем глобальные координаты из локальных)
	void Mapping(double *in, double *out);

	void Calc_J(int n_of_point); // вычисляет матрицу Якоби в точке Гаусса 
	void Calc_J(double x, double y, double z); // вычисляет матрицу Якоби в произвольной точке 
	void Calc_J_on_face(int n_of_point); // вычисляет Якобиан в точках, расположенных в грани
	void Calc_J_in_parallelepiped(); // вычисляет матрицу Якоби в случае параллелепипеда

	void Calc_local_matrix_b_for_hexahedron(); // матрица жёсткости (численно)
	void Calc_local_matrix_c_for_hexahedron(); // матрица массы (численно)

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

	// выдать x,y-компоненту поля в точках Гаусса в верхней грани шестигранника
	// сразу для трёх решений (для 3-слойной схемы по времени)
	void Get_x_y_on_face(double *ves1, double *ves2, double *ves3,
		double *outx1, double *outx2, double *outx3,
		double *outy1, double *outy2, double *outy3,
		double *outz1, double *outz2, double *outz3);

	// выдать z-компоненту ротора в точках Гаусса в верхней грани шестигранника
	// сразу для трёх решений (для 3-слойной схемы по времени)
	void Get_rotz_on_face(double *ves1, double *ves2, double *ves3,
		double *out1, double *out2, double *out3);

	// выдать компоненты ротора в точках Гаусса в верхней грани шестигранника
	// сразу для трёх решений (для 3-слойной схемы по времени)
	void Get_rot_on_face(double *ves1, double *ves2, double *ves3,
						  double *out1, double *out2, double *out3,
						  double *outx1, double *outx2, double *outx3,
						  double *outy1, double *outy2, double *outy3);

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

	void ScalarFieldOnParCff(double x, double y, double z, double *cff);

	//  Выдать в центре шестигранника, q - вектор весов
	double GetValueInHexCenter(double *q);
	void GetGradInHexCenter(double *q, double *out, int cj);

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


	// для ГЭЛ

	// выдать потенциал в точках Гаусса в верхней грани шестигранника
	// сразу для трёх решений (для 3-слойной схемы по времени)
	void GetAonFace(double *ves1, double *ves2, double *ves3,
		double *out1, double *out2, double *out3);

	// для гармонической петли
	void Set_dpr(double dpr);
	void Set_dpr0(double dpr0);
	void Set_mu0(double mu0);
	Vec_Prep_Data *d;


	//  для вычисления вектора правой части с учётом, что ток=const на элементе

	double asin0n[8][3], acos0n[8][3]; // нормальное поле в узлах
	double asin0c[3], acos0c[3]; // нормальное поле в центре элемента

	// вычислить локальный вектор правой части с учётом, что ток=const на элементе
	void LocalVectHarmConst();

	// вычислить нормальное поле в вершинах (для гармонических задач)
	void Calc_asin_acos_at_nodes();

	void GetVectorFieldNodes(double *ves, double *ax, double *ay, double *az);

	int npls,ipls,n_edges;

	void Calc_ss0();
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