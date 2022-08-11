/**                                                                                           
 * GENERAL REMARKS                                                                            
 *                                                                                            
 *  This code is freely available under the following conditions:                             
 *                                                                                            
 *  1) The code is to be used only for non-commercial purposes.                               
 *  2) No changes and modifications to the code without prior permission of the developer.    
 *  3) No forwarding the code to a third party without prior permission of the developer.     
 *                                                                                            
 *  			MTCalc_with_DFP_COCR                                                  
 *  Building a spline on the ground surface to output the value of the required function at an arbitrary point (receiver position)                                                                
 *                                                                                                                                                                                        
 *  Written by Ph.D. Petr A. Domnikov                                                         
 *  Novosibirsk State Technical University,                                                   
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                         
 *  p_domnikov@mail.ru                                                                        
 *  Version 1.2 March 10, 2021                                                                
*/                                                                                            

#pragma once
#include "AbstractFEM.h"
#include "OutputResultant3d.h"

class Give_out_vec_mt : public AbstractFEM3D
{
public:
	Vec_Prep_Data *d;

	T_Mapping_Vec *tmap;

	double *v3dat; // solution


	OutputResultant3d *resultantA,*resultantB;

	
	// normal field 
	std::complex <double> (*H_1d)[2];
	std::complex <double> (*E_1d)[2];
	// anomalous field in each of the receivers
	std::complex <double> (*H)[3];
	std::complex <double> (*E)[3];
	// impedance module in each of the receivers
	double *impedance;
	// apparent resistivity in each of the receivers
	double *rho;


	Give_out_vec_mt(Vec_Prep_Data *d, T_Mapping_Vec *tmap, double *v3dat);
	~Give_out_vec_mt();

	void Compute_1d_field();

	void Write_result_to_files();
	void Write_result_to_files_for_harm_loop();
	void Write_B_to_files_for_harm_loop(int StartType,int fdirect);
	void Write_E_to_files_for_harm_loop(int StartType,int fdirect);

	
	void Give_out_on_hex(bool for_harm_loop=false);

	int Read_1d_field(char *fname_in, double &Ex_s, double &Ex_c, double &Ey_s, double &Ey_c, 
		double &Hx_s, double &Hx_c,  double &Hy_s, double &Hy_c);
	double dA_dt(double t, double u_j, double u_j1, double u_j2,
		double dt, double dt0, double dt1, double t_j, double t_j1, double t_j2);


	int GetNumberOfNodes();
	int GetNumberOfElements();
	int GetElementNodesNumber();
	const pv::Point3D GetNode(const int& i_node);
	const pv::Point3D GetNodeTrue(const int& i_node);
	int GetNodeNumberOnElement(const int& i_element, const int& i_node);
	int GetElementMaterial(const int& i_element);
	int GetTypeOfElement(const int& i_element);
	double GetValueInElementCenter(const int& i_element, const Res3DValueType& r_type);
	int GetNumberOfResPoints(const Res3DValueType& r_type);
	pv::Point3D GetResPoint(const Res3DValueType& r_type, const int& i_point);
	int * GetPointerToRegular();
	int GetXSize();
	int GetYSize();
	int GetZSize();
	double *GetPointerToX();
	double *GetPointerToY();
	double *GetPointerToZ();
	void SaveResult(const Res3DValueType& r_type, const double& r_value, const int& i_point, const int& i_time);
};