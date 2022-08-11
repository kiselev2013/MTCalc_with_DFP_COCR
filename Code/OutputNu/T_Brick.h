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
 *  This file contains some basic routines for generation of local matrices and local vector of the right-hand side.     
 *  Calculation of values of basis functions inside a finite element.
 *  The T_Brick class contains functions for working with vector (edge-elements) for hexahedron mesh:
 * 	- calculation of local stiffness and mass matrices;
 * 	- calculation of the vector of the right-hand side through E_normal;
 *	- calculation of the Jacobi matrix at Gauss points and at an arbitrary point inside the hexahedron;
 * 	- output a solution and a curl inside a hexahedron or parallelepiped;
 * Nodal basis functions are needed here to calculate the transformation.
 * Template element [-1, 1]^3. (for both nodal and vector).
 * Vector basis functions with tangential components 2/hx, 2/hy, 2/hz along edges.                                                  
 *                                                                                                                       
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova,  Ph.D. Petr A. Domnikov, Ph.D. Yulia I. Koshkina      
 *  Novosibirsk State Technical University,                                                                              
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                                    
 *  Corresponding author:                                                                                                
 *  E-mail: p_domnikov@mail.ru (Petr A. Domnikov)                                                                        
 *  Version 1.3 March 30, 2021                                                                                           
*/                                                                                                                       

#pragma once
#include "ListForLine.h"
#include "Vec_Prep_Data.h"

//------------------------------------------------------------------------------------
class T_Brick
{
public:
	double hx, hy, hz; // element dimensions 
	double xk, xk1, yk, yk1, zk, zk1; //parallelepiped boundaries
    long num; // number of the final element in the grid

	double b[12][12]; // local stiffness matrix        
	double c[12][12]; // local mass matrix             
	double c0[12][12]; 	// local mass matrix
	double cSigma[12][12]; 	// local mass matrix	
	double cSigma0[12][12];	// local mass matrix 

	double f_re[12];  // values of the vector of the right-hand side in the midpoints of the edges

	double a[12][12]; // local matrix element    
	double g[12];     // local vector of element 
	double g_harm[24];
	double g_re[12], g_im[12];
	double g_re_sig[12], g_im_sig[12];
	double g_re_b[12], g_im_b[12];
	double g8[8]; // for nodes

	// equation coefficients
	double mu;
	double mu0;
	double sigma;
    Tensor sigmaTensor, sigmaTensor0, s_s0;
	double sigma0;
	double dpr;
	double dpr0;
	long n_mat; // material number of the element  

	long (*nver)[14];  // numbers of finite element nodes (hexahedrons with 13-nodes)   
	long (*ed)[25];    // elements listed by their edges + terminal edges + element type
	long *nvkat;       // material numbers of the finite elements                       
	long (*edges)[2];  // edges defined by 2 vertices                                   
	double (*xyz)[3];  // coordinates of the nodes                                      

	double *En; // normal field (on the edges)
	double (*En_nodes)[3]; // in nodes 

	//------------------start--AV formulation--start-------------------------------------------
	long (*nvetr)[20]; // numbers of nodes (1-8) and edges (9-20) in the nded array
	double MtrMass[20][20];
	double MtrMassRight[20][12];
	double MtrGest[20][20];
	double VctRight[40];
	double g_re_sig_av[20], g_im_sig_av[20];

	// constructor for MT
	T_Brick(long num, long (*nvetr)[20], long (*nver)[14], long (*ed)[25], long (*edges)[2], double (*xyz)[3], long *nvkat,
		Tensor *sigma3d, double *sigma0, double *mu3d, double omega,
		long alpha, long n_1d, double *z_1d, double *sin_1d, double *cos_1d);
	void Calc_block_local_matrix_and_vector_AV();
	// calculate the local stiffness matrix on the parallelepiped
	void Compute_Local_Matrix_B_AV(); 
	// calculate the local mass matrix on the parallelepiped  
	void Compute_Local_Matrix_C_AV(); 
	// compute a local 20x12 matrix to form the right side
	void Compute_Local_Matrix_Pr_AV(); 
	// calculate the local vector of the right side with the conversion of the 1-dimensional field into a 3-dimensional mesh for MT
	void Calc_local_vector_for_MT_AV();
	//------------------end--AV formulation--end-------------------------------------------

	int tasktype;

	// this constructor is usually used to return from an element
	T_Brick(double *x_coords, double *y_coords, double *z_coords, long type_of_hex);

	// constructor for a non-stationary task (for calculating local matrices)
	T_Brick(long num, long (*nver)[14], long (*ed)[25], long (*edges)[2], double (*xyz)[3], long *nvkat,
		double *sigma3d, double *sigma0, double *mu3d, double *En);

	// constructor for MT
	T_Brick(long num, long (*nver)[14], long (*ed)[25], long (*edges)[2], double (*xyz)[3], long *nvkat,
		double *sigma3d, double *sigma0, double *mu3d, double omega,
		long alpha, long n_1d, double *z_1d, double *sin_1d, double *cos_1d);

	T_Brick(long num, long (*nver)[14], double (*xyz)[3]);

	~T_Brick();

	// calculate local matrix for non-stationary problem
	// (mass matrix and vector or just add stiffness matrix)
	void Compute_Local_Matrix_And_Vector(const long what_compute); 

	// calculate the local stiffness matrix on the parallelepiped
	void Compute_Local_Matrix_B(); 

	// calculate the local mass matrix on the parallelepiped
	void Compute_Local_Matrix_C(); 

	// calculation of the vector of the right-hand side through E_Normal
	void Compute_Local_Vector_For_Anomal_Problem();

	// calculation of the vector of the right-hand side through the normal field, in the case when all the coefficients are discontinuous
	void ComputeLocalVectorMuEpsSigma(double *An, double *d2An);

	//	for MT	 

	double omega; 			// cyclic frequency for the harmonic task      
	double asin0[12], acos0[12]; 	// normal field in the midpoints of the edges  

	// for MT 1D problem  
	long alpha;
	long n_1d;
	double *z_1d;
	double *sin_1d;
	double *cos_1d;

	// calculate the local vector of the right-hand side with the conversion of the 1-dimensional field into a 3-dimensional mesh for MT
	void Calc_local_vector_for_MT();

	// calculate the normal field at the midpoints of the edges (for MT)
	void Calc_asin_acos_at_middle_of_edges();

	void Calc_block_local_matrix_and_vector();

	void Calc_J_Node_2(double x, double y, double z);
	double dPhi_node_2(long i, long j, double x, double y, double z);
	void Calc_V_Node_2(double *q,double x, double y, double z);
	double V[3];

	// for vector hexahedrons

	double x[8], y[8], z[8]; // coordinates of hexahedron vertices                   
	double J[3][3];          // Jacobian matrix                                      
	double J_1[3][3];	 	 // J^{-1}                                       
	double J_1_T[3][3];	 	 // J^{-T}                                       
	double det_J;            // Jacobian (determinant of the Jacobian matrix)        
	double det_J_abs;        // Jacobian modulus                                     
	long type_of_hex;        // 0..30 - parallelepiped, 31..61 - hexahedron
	double phi_all[12][3];    // basis functions at the point of integration              
	double rot_all[12][3];    // curls from basis functions at the point of integration   

	// converting a hexagon from template to arbitrary
	// (get global coordinates from local)
	void Mapping(double *in, double *out);

	void Calc_J(int n_of_point); 			// calculates the Jacobi matrix at the Gauss point                    
	void Calc_J(double x, double y, double z); 	// calculates the Jacobian matrix at an arbitrary point               
	void Calc_J_on_face(int n_of_point); 		// calculates the Jacobian at the points located in the face          
	void Calc_J_in_parallelepiped(); 		// calculates the Jacobian matrix in the case of parallelepiped       

	void Calc_local_matrix_b_for_hexahedron(); // stiffness matrix (numerically)
	void Calc_local_matrix_c_for_hexahedron(); // mass matrix (numerically)

	// Calculates the value of the vector field inside an arbitrary hexahedron
	void Calc_value_inside_hex(double *ves, double *in, double *out); 

	// value of the i-th basis function on the parallelepiped
	void Basis_func_on_vec_par(long i, double *in, double *out);

	// value of the i-th basis function on a parallelepiped with tangential components 2/hx, 2/hy, 2/hz
	void Basis_func_on_vec_par(long i, double ves, double *in, double *out);

	// value of the i-th basis function on the template parallelepiped 
	void Basis_func_on_reference_vec_par(long i, double *in, double *out);	

	// value of the i-th basis function on the hexahedron
	void Basis_func_on_vec_hex(long i, double ves, double *in, double *out);

	// output the curl at an arbitrary point inside the hexahedron
	void Calc_rotor_inside_hex(double *ves, double *in, double *out); 

	// the value of the curl of the i-th basis function on the template parallelepiped  
	void Rot_of_basis_func_on_reference_vec_par(long i, double *in, double *out);

	// the value of only the x-component of the curl of the i-th basis function inside the parallelepiped     
	void Rotx_of_basis_func_on_reference_vec_par(long i, double x, double *out);

	// the value of only the y-component of the curl of the i-th basis function inside the parallelepiped       
	void Roty_of_basis_func_on_reference_vec_par(long i, double y, double *out);

	// the value of only the z-component of the curl of the i-th basis function inside the parallelepiped
	void Rotz_of_basis_func_on_reference_vec_par(long i, double z, double *out);

	// the value of the curl  of the i-th basis function inside the parallelepiped  
	void Rot_of_basis_func_on_vec_par(long i, double ves, double *in, double *out);

	// give the x,y-component of the field at Gauss points in the upper face of the hexahedron
	// for three solutions at once (for a 3-layer scheme in time)                          	
	void Get_x_y_on_face(double *ves1, double *ves2, double *ves3,
		double *outx1, double *outx2, double *outx3,
		double *outy1, double *outy2, double *outy3,
		double *outz1, double *outz2, double *outz3);

	// give the z-component of the curl at Gauss points in the upper face of the hexahedron
	// for three solutions at once (for a 3-layer scheme in time)                          
	void Get_rotz_on_face(double *ves1, double *ves2, double *ves3,
		double *out1, double *out2, double *out3);

	// give the components of the curl at Gauss points in the upper face of the hexahedron
	// for three solutions at once (for a 3-layer scheme in time)                          
	void Get_rot_on_face(double *ves1, double *ves2, double *ves3,
						  double *out1, double *out2, double *out3,
						  double *outx1, double *outx2, double *outx3,
						  double *outy1, double *outy2, double *outy3);

	// Template one-dimensional basis functions on [-1, 1]^3 
	double l0(double x);
	double l1(double x);

	// Template nodal basis functions on [-1, 1]^3
	double Phi_node(long i, double x, double y, double z);

	// derivatives of nodal template basic functions
	double dPhi_node(long i, long j, double x, double y, double z);

	// print the value of a vector field on a parallelepiped 
	// (point coordinates - global)                          
	void VectorFieldOnPar(double x, double y, double z, double *ves,
		double *x_out, double *y_out, double *z_out);
	double ScalarFieldOnPar(double x, double y, double z, double *ves);
	double DxOfScalarFieldOnPar(double x, double y, double z, double *ves);
	double DyOfScalarFieldOnPar(double x, double y, double z, double *ves);
	double DzOfScalarFieldOnPar(double x, double y, double z, double *ves);

	void ScalarFieldOnParCff(double x, double y, double z, double *cff);

	// Output at the center of the hexahedron, q is the vector of weights
	double GetValueInHexCenter(double *q);
	void GetGradInHexCenter(double *q, double *out, int cj);

	// return the x-component of the field from the parallelepiped for three weight values at once
	void VectorFieldXOnPar3(double y, double z, double *ves_j2, double *ves_j1, double *ves_j,
		double *out_j2, double *out_j1, double *out_j);

	// return the y-component of the field from the parallelepiped for three weight values at once
	void VectorFieldYOnPar3(double x, double z, double *ves_j2, double *ves_j1, double *ves_j,
		double *out_j2, double *out_j1, double *out_j);

	void RotXOnPar(double x, double *ves, double *out, bool loc_c=false);
	void RotYOnPar(double y, double *ves, double *out, bool loc_c=false);
	void RotZOnPar(double z, double *ves, double *out);
	void RotZOnPar3(double z,
		double *ves1, double *ves2, double *ves3,
		double *out1, double *out2, double *out3);

	// variable substitution
	void Transformation_of_variables(double *in, double *out);
	void Transformation_of_variables(double *x, double *y, double *z);
	double Xi(double x);
	double Eta(double y);
	double Zeta(double z);


	

	// give the potential at Gauss points in the upper face of the hexahedron 
	// for three solutions at once (for a 3-layer time scheme) 
	void GetAonFace(double *ves1, double *ves2, double *ves3,
		double *out1, double *out2, double *out3);

	
	void Set_dpr(double dpr);
	void Set_dpr0(double dpr0);
	void Set_mu0(double mu0);
	Vec_Prep_Data *d;


	// to calculate the vector of the right-habd side, taking into account that the current=const on the element

	double asin0n[8][3], acos0n[8][3]; //normal field at nodes
	double asin0c[3], acos0c[3]; // normal field in the center of the element

	// calculate the local vector of the right-hand side, taking into account that current=const on the element
	void LocalVectHarmConst();

	// calculate the normal field at the vertices (for harmonic problems)
	void Calc_asin_acos_at_nodes();

	void GetVectorFieldNodes(double *ves, double *ax, double *ay, double *az);

	int npls,ipls,n_edges;

	void Calc_ss0();
};
//------------------------------------------------------------------------------
const double MIDDLE_OF_LOCAL_EDGE[12][3] = {
	 0.0, -1.0, -1.0,
	 0.0,  1.0, -1.0,
	 0.0, -1.0,  1.0,
	 0.0,  1.0,  1.0,

	-1.0,  0.0, -1.0,
	-1.0,  0.0,  1.0,
	 1.0,  0.0, -1.0,
	 1.0,  0.0,  1.0,

	-1.0, -1.0,  0.0,
	 1.0, -1.0,  0.0,
	-1.0,  1.0,  0.0,
	 1.0,  1.0,  0.0
};
//-----------------------------------------------------------------------
const double LOCAL_COORDS_OF_NODES[8][3] = 
{
	-1.0, -1.0, -1.0,
	 1.0, -1.0, -1.0,
	-1.0,  1.0, -1.0,
	 1.0,  1.0, -1.0,

	-1.0, -1.0,  1.0,
	 1.0, -1.0,  1.0,
	-1.0,  1.0,  1.0,
	 1.0,  1.0,  1.0
};
//-----------------------------------------------------------------------
const double TANGENT_VECTORS_ON_REFERENCE_CUBE[12][3] = {
	1.0, 0.0, 0.0,
	1.0, 0.0, 0.0,
	1.0, 0.0, 0.0,
	1.0, 0.0, 0.0,

	0.0, 1.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 1.0, 0.0,

	0.0, 0.0, 1.0,
	0.0, 0.0, 1.0,
	0.0, 0.0, 1.0,
	0.0, 0.0, 1.0
};
//-----------------------------------------------------------------------
const long REG_EDGES[12][2]={ 	// what non-terminal edges are on the element
	0,1, 2,3, 4,5, 6,7, 0,2, 4,6, 1,3, 5,7, 0,4, 1,5, 2,6, 3,7 };
//-----------------------------------------------------------------------