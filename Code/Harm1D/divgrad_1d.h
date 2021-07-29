#pragma once
class Divgrad_1d
{
public:
	double *coords_1d; 
	double *sin_1d;
	double *cos_1d;
	long n_1d; // число узлов в одномерной сетке
	double bak; // нижн€€ граница бака
	double step0; // начальный шаг (у поверхности «емли) дл€ одномерной сетки
	double coef_razr; // коэффициент разр€дки
	long n_layers_1d; // число слоЄв
    double *layers_1d; // границы слоЄв
	double *sigma_1d;
	long alpha;

	int Solve_1d_Problem_for_3d_task();
	
	Divgrad_1d();
	~Divgrad_1d();
};
