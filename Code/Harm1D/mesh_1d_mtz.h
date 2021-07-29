#pragma once
const double X0_1D = -10000000.0; // нижняя граница бака
const double X1_1D = 0.0;
const double COEF_RAZR = 1.02; // разрядка
const double STEP0 = 1e-2;     // начальный шаг
//------------------------------------------------
class Mesh_1d_mtz
{
public:
	long n_points; 
	long n_elem; 
	long n_materials;
	long n_layers;
	double *coords;
	long *nvkat;
	double *layers;
	long *nodes_in_layer;

	double *mu;
	double *sigma;
	double omega;

	Mesh_1d_mtz();
	~Mesh_1d_mtz();

	int Read_1d_data_for_3d_problem();
	void Gen_1d_mesh();
};
