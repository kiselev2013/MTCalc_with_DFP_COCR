#pragma once

class Vec_Prep_Data
{
public:

	Vec_Prep_Data();
	~Vec_Prep_Data();

	int Read_prep_data(); // чтение сетки для МТЗ на векторных
	int Read_mtz_1d(); // считывает usin.dat, ucos.dat, alfa, nu (это для МТЗ)

	int Read_3dmeshregular(long interval); 

	int LoadReceivers(char *fname, int& _n_pointres, double (*(&_pointres))[3]);

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
	long n_layers_1d;  // кол-во различных слоёв (материалов) для одномерки (из sreda1d.ay)
	double *layers_1d; // слои для одномерки (из sreda1d.ay)
	double *sigma_1d;  // sigma0 для одномерки (из sreda1d.ay)

	int N_X,N_Y,N_Z;
	vector<int> regular;
	vector<double> Xcrd,Ycrd,Zcrd;

	long n_mesh_regular_x;  // число шагов по x в 3dmeshregular
	long n_mesh_regular_y;  // число шагов по y в 3dmeshregular
	double *mesh_regular_x; // x-координаты из 3dmeshregular
	double *mesh_regular_y; // y-координаты из 3dmeshregular

 	double nu;    // частота
 	long alfa;    // направление тока: по x - alpha=1; по y - alpha=0; J=(alpha*Jx,(1-alpha)*Jy, 0)
 	double *usin; // решение одномерной задачи (sin-компонента)
 	double *ucos; // решение одномерной задачи (sin-компонента)
 	double *z_1d; // одномерная сетка
 	long n_1d;    //  число узлов в одномерной сетке

	int npr, nfreq;
};
//-----------------------------------------------------------