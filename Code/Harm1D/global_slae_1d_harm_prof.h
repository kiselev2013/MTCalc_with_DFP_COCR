#pragma once
class Global_slae_1d_harm_prof
{
public:

	long n;       // размерность СЛАУ
	long n_elem;  // число элементов в сетке

	long *ig;    // массив адресов начал строк
	double *ggl; // элементы нижнего треугольника
	double *ggu; // элементы верхнего треугольника
	double *di;  // диагональ
	double *pr;  // вектор правой части

	long ig_n_1; // размерность ggl, ggu

	long *nvkat; // номера материалов

	double *xyz; // сетка

	long ay_hy;

	// коэффициенты уравнения
	double *mu;
	double *sigma;
	double omega;

	Global_slae_1d_harm_prof(long n, long n_elem, long *ig, long *nvkat,
		double *mu, double *sigma, double omega, double *xyz, long ay_hy);
	~Global_slae_1d_harm_prof();

	void Assembling_for_1d_harm_problem();
	void Set_boundary_conditions_for_Ay();
	void Set_boundary_conditions_for_Hy();
};
