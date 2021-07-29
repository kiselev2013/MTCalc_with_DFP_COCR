#pragma once
#include "T_Mapping.h"
//---------------------------------------------------------------------------------
/*
Класс Vec_Prep_Data служит для хранения сетки и другой информации о задаче,
которую GeoPrep выкладывает на диск в файлы.

Текущая реализация предназаначена для петли и МТЗ на векторных шестигранниках.
МТЗ на векторных пока не встроено, но функции для чтения сетки я решил оставить.

из файла 3dmeshregular пока считываются только x,y-координаты. (z не считывается)
*/
//---------------------------------------------------------------------------------
struct Tensor
{
	double val[3][3];

	Tensor()
	{
		Clear();
	}

	void Clear()
	{
		for(int i=0; i<3; i++)
		{
			for(int j=0; j<3; j++)
				val[i][j] = 0;
		}
	}

	Tensor& operator = (const double& d)
	{
		Clear();

		for(int j=0; j<3; j++)
			val[j][j] = d;

		return *this;
	}

	bool operator == (const double& d)
	{
		bool flag=true;
		int i, j;
		const double eps = 1e-6;

		for (i=0; i<3; i++)
		{
			for (j=0; j<3; j++)
			{
				if (i==j)
				{
					if(fabs(val[i][i] - d) > eps)
					{
						flag = false;
						break;
					}
				}
				else
				{
					if (val[i][j] != 0)
					{
						flag = false;
						break;
					}
				}
			}

			if (!flag)
				break;
		}

		return flag;
	}

	bool operator != (const double& d)
	{
		return !(*this == d);
	}
};
//------------------------------------------------------------------------

class Vec_Prep_Data
{
public:

	struct EnLine
	{
		long mtr;
		double es[3];
		double ec[3];
	};

	Vec_Prep_Data();
	~Vec_Prep_Data();

	int Read_prep_data(); // чтение сетки для МТЗ на векторных
	int ReadPrepDataHarmLoop(char *pointres_fname); // чтение сетки для петли с гармоническим источником
	int Read_mtz_1d(); // считывает usin.dat, ucos.dat, alfa, nu (это для МТЗ)

	int Read_mesh_for_nonstat_problem(char *pointres_fname); // для петли на векторных
	int Read_infite0();  // считывает кол-во времён и времена из infite.0

	// читает 3dmeshregular, пропуская координаты сетки с интервалом=interval
	// (если хотим построить сплайн с шагом больше, чем в 3dmeshregular)
	int Read_3dmeshregular(long interval); 

	int LoadReceivers(char *fname, int& _n_pointres, double (*(&_pointres))[3]);

	long n_pointres;
	double (*pointres)[3]; // координаты приёмников

	long maxiter; // максимальное число итераций для решателя (из config)
	double eps;   // epsilon для решателя (из config)
	long n_materials; // число различных материалов
	double *mu3d;     //(mu для трёхмерной задачи с объектами (первый столбик из файла mu3d))
	double *mu0;      //(mu для одномерки (второй столбик из файла mu3d))
	int n_pointresB;  // число приёмников для B
	int n_pointresE;  // число приёмников для E
	double (*pointresB)[3]; // координаты приёмников
	double (*pointresE)[3]; // координаты приёмников
	double *sigma3d;       //(sigma для трёхмерной задачи с объектами (первый столбик из sig3d))
	double *sigma0;        //(sigma для одномерки (второй столбик из sig3d))
	double *dpr3d;       //(dpr для трёхмерной задачи с объектами (первый столбик из dpr3d))
	double *dpr0;        //(dpr для одномерки (второй столбик из dpr3d))
	long kuzlov;      // число узлов (всех, в том числе и терминальных)
	long kpar;        // число элементов в сетке
	long kt1;         // число узлов с первыми краевыми
	long *l13d;       // номера узлов с первыми краевыми
	long (*nver)[14]; // номера узлов конечных элементов (13-узловые шестигранники)
	long *nvkat;       // номера материалов конечных эл-тов
	double (*xyz)[3];  // координаты узлов
	double (*xyzt)[3];  // координаты узлов
	long n_layers_1d;  // кол-во различных слоёв (материалов) для одномерки (из sreda1d.ay)
	double *layers_1d; // слои для одномерки (из sreda1d.ay)
	double *sigma_1d;  // sigma0 для одномерки (из sreda1d.ay)

	int LoadVectorE0ForLine(int n);
	vector< vector<EnLine> > EnForLine;

	int N_X,N_Y,N_Z;
	vector<int> regular;
	vector<double> Xcrd,Ycrd,Zcrd;

	// 3dmeshregular
	long n_mesh_regular_x;  // число шагов по x в 3dmeshregular
	long n_mesh_regular_y;  // число шагов по y в 3dmeshregular
	double *mesh_regular_x; // x-координаты из 3dmeshregular
	double *mesh_regular_y; // y-координаты из 3dmeshregular

	// только для нестационарной задачи
	long ntime;     // всего временных слоёв
	double *time;   // временные слои

	// только для МТЗ 
 	//long norvect; // число векторов для GMRES
 	double nu;    // частота
 	long alfa;    // направление тока: по x - alpha=1; по y - alpha=0; J=(alpha*Jx,(1-alpha)*Jy, 0)
 	double *usin; // решение одномерной задачи (sin-компонента)
 	double *ucos; // решение одномерной задачи (sin-компонента)
 	double *z_1d; // одномерная сетка
 	long n_1d;    //  число узлов в одномерной сетке

	T_Mapping_Vec *tmap;
	long tasktype;

	int npr, nfreq;

	int fdirect;
};
//-----------------------------------------------------------