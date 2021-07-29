#pragma once
//------------------------------------------------------------------------
	double Scal(double *a, double *b, long n);
	double Norm_Euclid(double *a, long n);
	double Norm_Max(double *a, long n);
	double Projection_On_Axis(double *v,double *o); //Проекция вектора v на ось o
	void Mult_Plot(double *a, double *x, double *y, long n);
	void Mult_MV(long *ig, long *jg, double *ggl, double *ggu, double *di, double *x, double *y, long n);
	double Relative_Error(double *analytic, double *numeric, long n);
	long Max_Long(long a, long b);
	long Min_Long(long a, long b);
	void Sort2(long *a, long *b);
	double Interval(double *x, double *y);
	double Interval_Parallel_Lines(double *a0, double *a1, double *b0, double *b1);
	double Spline(double x, long n, double *xyz, double *values);

	double Calc_dof(double *J, double *func, long n_local_edge);

	// сообщает об ошибке выделения памяти и генерирует исключение
	void Memory_allocation_error(const char *var, const char *func);

	// сообщает об ошибке открытия файла и генерирует исключение
	void Cannot_open_file(const char *fname, const char *func);

	// сообщает об ошибке открытия файла, но работа программы продолжается
	void Cannot_open_file_but_continue(const char *fname, const char *func);

	void Mult_Plot_AV(double *a, double *x, double *y, long n, long m);

