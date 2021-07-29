#pragma once
//--------------------------------------------
class Hex_Local_Matrix
{
public:
	int number_of_element;
	double x[8], y[8], z[8]; // координаты вершин шестигранника	
	double J[3][3];          // матрица Якоби
	double J_1_T[3][3];		 // J^{-T}
	double det_J;            // Якобиан (определитель матрицы Якоби)
	double det_J_abs;        // модуль Якобиана
	int type_of_hex; // 0..30 - параллелепипед, 31..61 - шестигранник
	double grad_all[8][3]; // градиенты от баз. функций в точке интегрирования
	// локальные матрицы 8х8 (обычные, не блочные)
	double b[8][8]; // матрица массы
	double hx, hy, hz; // для параллелепипеда
	bool JforParCalc; // вычислялся ли ранее Якобиан в случае параллелепипеда
	void Calc_J(int n_of_point);
	void Calc_local_matrix_b_for_parallelepiped();
	void CalcMassMatrix();
	Hex_Local_Matrix(int i, long (*nver)[14], double (*xyz)[3]);
	~Hex_Local_Matrix();
};
