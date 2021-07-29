#pragma once
//------------------------------------------------------------------------
	double Scal(double *a, double *b, int n);
	double Norm_Euclid(double *a, int n);
	double Norm_Max(double *a, int n);
	double Projection_On_Axis(double *v,double *o); //Проекция вектора v на ось o
	void Mult_Plot(double *a, double *x, double *y, int n);
	void Mult_MV(int *ig, int *jg, double *ggl, double *ggu, double *di, double *x, double *y, int n);
	double Relative_Error(double *analytic, double *numeric, int n);
	int Max_Long(int a, int b);
	int Min_Long(int a, int b);
	void Sort2(int *a, int *b);
	double Interval(double *x, double *y);
	double Interval_Parallel_Lines(double *a0, double *a1, double *b0, double *b1);
	double Spline(double x, int n, double *xyz, double *values);

	// сообщает об ошибке выделения памяти и генерирует исключение
	void Memory_allocation_error(const char *var, const char *func);

	// сообщает об ошибке открытия файла и генерирует исключение
	void Cannot_open_file(const char *fname, const char *func);

	// сообщает об ошибке открытия файла, но работа программы продолжается
	void Cannot_open_file_but_continue(const char *fname, const char *func);

	// для выделения sin-, cos-компонент из нестационарной задачи в одной точке
	int Spline_sin_cos(double *x, double *t, int m, double w, double *e);

//	int Read_EBxyz(char *fname, Mesh3DWithNormalField* ptrMesh3D, int beg, int end);
//	int Read_EBxyz_all(Mesh3DWithNormalField* ptrMesh3D);
	void Mult_Plot_AV(double *a, double *x, double *y, int n, int m);

