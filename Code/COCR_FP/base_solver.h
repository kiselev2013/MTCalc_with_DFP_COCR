#pragma once
//------------------------------------------------------------------------
class Base_solver
{
//  ласс Base_solver содержит элементарные подпрограммы дл€ построени€ 
// решателей (скал€рное произведение, норма вектора и т.д.)
public:
	Base_solver();
	~Base_solver();

	inline double Scal(double *x, double *y, int n);     // скал€рное произведение
	double Norm_Euclid(double *x, int n);         // евклидова норма вектора
	double Projection(double *vec, double *axis);  // проекци€ вектора на ось
	double Relative_Error(double *analytic, double *numeric, int n); // относительна€ погрешность
	double Spline(double x, int n, double *xyz, double *values);

	//”множение матрицы (в плотном формате) на вектор
	void Mult_Plot(double *a,double *pr,double *rez,int n); 

	// вращени€ √ивенса
	int  Givens1(double& x, double& y, double& c, double& s);
	void Givens2(double& x, double& y, double c, double s);
	int  Givens(double *a, double *f, int n);

	// обратный ход по плотной матрице
	int Undirect(double *a, double *b, double *x, int n);

	// –ешение —Ћј” с квадратной матрицей, нижний треугольник которой
	// содержит только одну ненулевую поддиагональ.
	// Ёто необходимо, если вылетит ортогонализаци€ јрнольди
	int Solve_square_subdiag(double *a, double *b, double *x, int n);

	// запись в файл относительной нев€зки, с которой вышли, eps,
	// числа итераций и времени решени€ —Ћј”
	int WriteKitChrono(char *fname, double residual, double eps, int iter, double time);
	int Write_kit(char *fname, double residual, double eps, int iter, int time);
	int Write_kit(char *fname, double residual, double eps, int iter, int time, int n);
	int Write_kit(char *fname, double residual, double eps, int iter, int time, int n, int size_jg);
	int Write_kit(char *fname, double residual, double eps, int iter, int time, double change_of_solution);


	// установить число итераций дл€ решени€ —Ћј” с факторизованнной (под)матрицей
	void SetPrcndItr(int n);

	// вернуть относительную норму нев€зки с последней итерации
	double GetResidual();

	// запомнить относительную норму нев€зки с последней итерации
	void SetResidual(double residual);



	// дл€ комлекснозначной арифметики (векторы хран€тс€ в обычном формате, комплексные числа (скал€ры) - в std::complex<double>)

	// скал€рное произведение дл€ комплекснозначных векторов
	std::complex<double> ScalCmplxTrue(double *x, double *y, int nb); 

	// комплексно-сопр€женное скал€рное произведение
	std::complex<double> ScalCmplx(double *x, double *y, int nb); 

	// умножение вектора на комплексное число
	void MultCmplxNumVect(std::complex<double> a, double *x, double *y, int nb);

	// умножение компонент одного комплексного вектора на компоненты другого комплексного вектора
	void MultCmplxVectVect(int nb, double *a, double *b, double *c);

	// деление компонент одного комплексного вектора на компоненты другого комплексного вектора
	void DivCmplxVectVect(int nb, double *a, double *b, double *c);

	// z = x + a*y
	void Cmplx_axpy(std::complex<double> a, double *x, double *y, double *z, int nb);

	double GetMinimiz(int n, double *ap, double *r);
};

//----------------------------------------------------
