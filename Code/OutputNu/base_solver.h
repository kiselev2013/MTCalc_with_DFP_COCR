#pragma once

class Base_solver
{
/*
 ласс Base_solver содержит элементарные подпрограммы дл€ построени€ 
решателей (скал€рное произведение, норма вектора и т.д.)
*/

// итерации в решател€х нумеруютс€ с 1

protected:
	double *y_omp;


public:
	Base_solver();
	~Base_solver();

	int n_threads; // число потоков (процессоров дл€ OpenMP)

	double Scal(double *x, double *y, long n);     // скал€рное произведение
	double Norm_Euclid(double *x, long n);         // евклидова норма вектора
	double Projection(double *vec, double *axis);  // проекци€ вектора на ось
	double Relative_Error(double *analytic, double *numeric, long n); // относительна€ погрешность
	double Spline(double x, long n, double *xyz, double *values);

	//”множение матрицы (в плотном формате) на вектор
	void Mult_Plot(double *a,double *pr,double *rez,long n); 

	// вращени€ √ивенса
	int  Givens1(double& x, double& y, double& c, double& s);
	void Givens2(double& x, double& y, double c, double s);
	int  Givens(double *a, double *f, long n);

	// обратный ход по плотной матрице
	int Undirect(double *a, double *b, double *x, long n);

	// –ешение —Ћј” с квадратной матрицей, нижний треугольник которой
	// содержит только одну ненулевую поддиагональ.
	// Ёто необходимо, если вылетит ортогонализаци€ јрнольди
	int Solve_square_subdiag(double *a, double *b, double *x, long n);

	// запись в файл относительной нев€зки, с которой вышли, eps,
	// числа итераций и времени решени€ —Ћј”
	int Write_kit(char *fname, double residual, double eps, long iter, __time64_t time);
	int Write_kit(char *fname, double residual, double eps, long iter, __time64_t time, long n);
	int Write_kit(char *fname, double residual, double eps, long iter, __time64_t time, long n, long size_jg);
	int Write_kit(char *fname, double residual, double eps, long iter, __time64_t time, double change_of_solution);


	// дл€ комлекснозначной арифметики (векторы хран€тс€ в обычном формате, комплексные числа (скал€ры) - в std::complex<double>)

	// скал€рное произведение дл€ комплекснозначных векторов
	std::complex<double> ScalCmplx(double *x, double *y, long nb); 

	// умножение вектора на комплексное число
	void MultCmplxNumVect(std::complex<double> a, double *x, double *y, long nb);

	// умножение компонент одного комплексного вектора на компоненты другого комплексного вектора
	void MultCmplxVectVect(long nb, double *a, double *b, double *c);

	// деление компонент одного комплексного вектора на компоненты другого комплексного вектора
	void DivCmplxVectVect(long nb, double *a, double *b, double *c);

protected:
	long nMRSrestart; // число итераций, через которое перезапускаетс€ MRS
	// (чтобы не накапливались ошибки округлени€)
public:
	void Set_nMRSrestart(long nMRSrestart) {this->nMRSrestart = nMRSrestart;}

};
//----------------------------------------------------
