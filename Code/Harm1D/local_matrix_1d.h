#pragma once
class Local_matrix_1d
{
public:

	double b[2][2]; // матрица жёсткости
	double c[2][2]; // матрица массы 
	double f_re[2], f_im[2]; // значения правой части в узлах
	double g_re[2], g_im[2]; // векторы правой части для Re и Im уравнений
	double a[4][4]; // матрица для гармонической задачи
	double g[4]; // вектор правой части для гармонической задачи

	// коэффициенты уравнения
	double mu;
	double sigma;
	double omega;

	double h; // размер конечного элемента

	long alpha; // направление тока


	Local_matrix_1d(double h, double mu, double sigma, double omega,
		double *f_re, double *f_im, long alpha);
	~Local_matrix_1d();

	void Compute_local_matrix_and_vector_harm();

	void Compute_local_matrix_harm();
	void Compute_local_vector_harm();
	void Compute_local_matrix_b();
	void Compute_local_matrix_c();
	void Compute_local_vector(double *f, double *g);
};
