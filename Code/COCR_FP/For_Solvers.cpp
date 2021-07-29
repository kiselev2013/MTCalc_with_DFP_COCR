#include "stdafx.h"
#ifdef GEOPREP_MPI
#include "..\Parallel\ClusterMPI.h"
	extern int rankMPI;
	extern int sizeMPI;
	extern int rankClient;
#endif
extern ofstream logfile;
//------------------------------------------------------------------------
// сообщает об ошибке выделения памяти
//------------------------------------------------------------------------
void Memory_allocation_error(const char *var, const char *func)
{
	string str;
	str = "MEMORY ALLOCATION ERROR for variable ";
	str = str + '\"' + var + '\"' + " in function " + '\"' + func + "\"\n";
	logfile << str << flush;
	cerr    << str << flush;
	cout << str << flush;
	throw logic_error(str);
}
//------------------------------------------------------------------------
// сообщает об ошибке открытия файла и генерирует исключение
//------------------------------------------------------------------------
void Cannot_open_file(const char *fname, const char *func)
{
	string str;
	str = "CANNOT OPEN FILE ";
#ifdef GEOPREP_MPI
	if (GEOPREP_MPI_SERVER)
	{
		str += "(Server) ";
	} 
	else
	{
		str += "(Client) ";
	}
#endif
	str = str + '\"' + fname + '\"' + " in function " + '\"' + func + "\"\n";
	logfile << str << flush;
	cerr    << str << flush;
	cout << str << flush;
	throw logic_error(str);
}
//------------------------------------------------------------------------
// сообщает об ошибке открытия файла, но работа программы продолжается
//------------------------------------------------------------------------
void Cannot_open_file_but_continue(const char *fname, const char *func)
{
	string str;
	str = "Cannot open file ";
	str = str + '\"' + fname + '\"' + " in function " + '\"' + func + "\"\n";
	logfile << str << flush;
	cerr    << str << flush;
}
//------------------------------------------------------------------------
// для выделения sin-, cos-компонент из нестационарной задачи в одной точке
//------------------------------------------------------------------------
double Scal(double *a, double *b, int n)
{
	int i;
	double sum = 0.0;

	for(i=0; i<n; i++)
		sum += a[i]*b[i];

	return sum;
}
//------------------------------------------------------------
void Mult_Plot_AV(double *a, double *x, double *y, int n, int m)
{
	int i, j, temp;
	double sum;

	for(i=0;i<n;i++)
	{
		sum = 0.0;
		temp = i*m;
		for(j=0;j<m;j++)
			sum += a[temp+j]*x[j];

		y[i] = sum;
	}
}
//------------------------------------------------------------
double Projection_On_Axis(double *v,double *o) //Проекция вектора v на ось o
{
	double value;

	value = Scal(v,o,3)/sqrt(Scal(o,o,3));

	return value;
}
//------------------------------------------------------------
void Mult_Plot(double *a, double *x, double *y, int n)
{
	int i, j, temp;
	double sum;

	for(i=0;i<n;i++)
	{
		sum = 0.0;
		temp = i*n;
		for(j=0;j<n;j++)
			sum += a[temp+j]*x[j];

		y[i] = sum;
	}
}
//------------------------------------------------------------
double Norm_Euclid(double *a, int n)
{
	double value;

	value = sqrt(Scal(a,a,n));

	return value;	
}
//------------------------------------------------------------
double Norm_Max(double *a, int n)
{
	double max, current;
	int i;

	max = 0.0;

	for(i=0; i<n; i++)
	{
		current = a[i]*a[i];
		if(current > max)
			max = current;
	}

	return max;
}
//-------------------------------------------------------------
//---------       умножение матрицы на вектор        ----------
//-------------------------------------------------------------
void Mult_MV(int *ig, int *jg, double *ggl, double *ggu, double *di, double *x, double *y, int n)
{
	int i, j, k;

	for(i=0; i<n; i++)
	{
		y[i] = di[i]*x[i];
		for(j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			y[i] += ggl[j]*x[k];
			y[k] += ggu[j]*x[i];
		}
	}
}
//-----------------------------------------------------------
double Relative_Error(double *analytic, double *numeric, int n)
{
	double *razn=NULL;
	double norm_analytic;
	double norm_razn;
	double error;
	int i;

	razn = new double[n];
	if(razn == 0)
		Memory_allocation_error("razn", "Relative_Error");

	for(i=0;i<n;i++)
		razn[i] = analytic[i] - numeric[i];

	In_Out R;
	R.Write_Txt_File_Of_Double("razn", razn, n, 1);

	norm_analytic = Norm_Euclid(analytic, n);
	norm_razn = Norm_Euclid(razn, n);

	error = norm_razn/norm_analytic;

	if(razn) {delete [] razn; razn=NULL;}

	return error;
}
//-----------------------------------------------------------
int Max_Long(int a, int b)
{
	if(a > b) return a;
	return b;
}
//-----------------------------------------------------------
int Min_Long(int a, int b)
{
	if(a < b) return a;
	return b;
}
//-----------------------------------------------------------
void Sort2(int *a, int *b)
{
	int tmp;
	if(*a < *b)
	{
		tmp = *a;
		*a = *b;
		*b = tmp;
	}
}
//-----------------------------------------------------------
double Interval(double *x, double *y)
{
	double r;

	r = sqrt((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]) + (x[2]-y[2])*(x[2]-y[2]));

	return r;
}
//-----------------------------------------------------------
double Interval_Parallel_Lines(double *a0, double *a1, double *b0, double *b1)
{
	double ax, ay, az; // направляющий вектор для 1-й прямой
	double x1, y1, z1; // точка на первой прямой
	double x0, y0, z0; // точка на второй прямой
	double det1, det2, det3;
	double answer;

	ax = a1[0] - a0[0];
	ay = a1[1] - a0[1];
	az = a1[2] - a0[2];

	x1 = a0[0];
	y1 = a0[1];
	z1 = a0[2];

	x0 = b0[0];
	y0 = b0[1];
	z0 = b0[2];

	det1 = ay*(z1 - z0) - az*(y1 - y0);
	det2 = az*(x1 - x0) - ax*(z1 - z0);
	det3 = ax*(y1 - y0) - ay*(x1 - x0);

	answer = sqrt(det1*det1 + det2*det2 + det3*det3)/sqrt(ax*ax + ay*ay + az*az);

	return answer;
}
//-----------------------------------------------------------
double Spline(double x, int n, double *xyz, double *values)
{
	// n - число элементов
	double s, xi;
	int i, t, flag;

	flag = 0;

	if (x < xyz[0])
	{
		return values[0];
	}

	if (x > xyz[n-1])
	{
		return values[n-1];
	}

	for(i=0; i<n-1; i++)
	{
		if(x >= xyz[i]  &&  x <= xyz[i+1])
		{
			t = i;
			flag = 1;
			break;
		}
	}

	if(flag == 1)
	{
		xi = (x - xyz[t])/(xyz[t+1] - xyz[t]);
		s = (1.0 - xi)*values[t] + xi*values[t+1];
	}
	else
	{
		s = 0.0;
	}

	return s;
}
