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
 *  This file contains implementation of the class for transition matrix used to store a nonconforming mesh
 *                                                                                                         
 *                                                                                                         
 *  Written by Ph.D. Petr A. Domnikov                                                                      
 *  Novosibirsk State Technical University,                                                                
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                      
 *  p_domnikov@mail.ru                                                                                     
 *  Version 2.1 March 17, 2021                                                                             
*/                                                                                                         


#pragma once

//-------------------------------------------------------------------------
// The long_double structure is required when constructing the T-matrix.
// It stores the number of the basic function and the value with which it is taken.
// It is convenient to store these two numbers together.
//-------------------------------------------------------------------------
struct long_double
{
	long i;
	double d;
};
//-------------------------------------------------------------------------
// A binary tree is used to store "chains" of numbers of basis functions when constructing a T-matrix 
//-------------------------------------------------------------------------
class Btree
{
public:
	long elem_long;
	double elem_double;
	bool visited;
	Btree *left, *right;
	
	Btree();
	Btree(long elem_long, double elem_double);
	~Btree();

	int DeleteTree(Btree *t);
	int Add_Left(long elem_long, double elem_double);
	int Add_Right(long elem_long, double elem_double);
	Btree * visit(std::vector<long_double> *s_t, Btree *t, long j);
};
//--------------------------------------------------------------------------------------------------------------
// The T_Mapping_Vec class is used to construct a T-matrix for vector (edge-elements) hexagons and numbering edges in a nodal grid.                                                                                                           
// Here we use vector basis functions with tangential components 2/hx, 2/hy, 2/hz along the edges,
// therefore, the sum over the column in the T-matrix is not equal to 1, and such a check is meaningless. 
//--------------------------------------------------------------------------------------------------------------
class T_Mapping_Vec
{
public: 
	T_Mapping_Vec(long (*nver)[14], double (*xyz)[3], long kuzlov, long kpar);
	~T_Mapping_Vec();

	// transition matrix in sparse column format
	long *ig_t; 
	long *jg_t;
	double *gg_t;

	long *ig_s, *jg_s; // SIGMA structure (auxiliary structure for T-matrix)
	double *s_val;     // values that are required to form all nonzero components of the matrix T 

	long n_c;  // number of continuous functions (continuous)
	long n_dc; // number of discontinuous functions (discontinuous)
	long n;    // the total number of functions is n = n_c + n_dc

	long kuzlov; // number of nodes in the grid
	long kpar;   // number of parallelepipeds     +
	long (*nver)[14]; // cells listed by their nodes + terminal nodes + element type
	double (*xyz)[3]; // vertex coordinates (all together)

	long (*edges)[2]; // edges defined by 2 vertices
	long (*ed)[25];   // cells listed by their edges + terminal edges + element type

	//--------------------------------------------------------------------------------------
	// for A-V formulation
	// for edge functions
	long edges_c;  // number of continuous functions (continuous)
	long edges_dc; // number of discontinuous functions (discontinuous)
	long edges_all;    // the total number of functions is edges_all = edges_c + edges_dc
	// for nodal functions
	long nodes_c;  // number of continuous functions (continuous)      
	long nodes_dc; // number of discontinuous functions (discontinuous)
	long nodes_all;    // the total number of functions is  nodes_all = nodes_c + nodes_dc
	// unknown
	long unk_all;

	long *nvkat;
	long *m_nded;		// m_nded - array for storing numbers of unknowns (nodes and edges)
	long *m_nded_type;	// 0 - node, 1 - edge
	long (*nvetr)[20];	// elements listed by their nodes and edges in the numbering corresponding to m_nded
	// 1-8  - node numbers in the array m_nded
	// 9-20 - edge numbers in m_nded array
	long *nodes_position_in_nded;	// positions of nodes in the m_nded array
	long *edges_position_in_nded;	// edge positions in m_nded array
	bool *nodes_earth; // 0 - the node is in the air, 1 - in the ground
	int Form_data_A_V(long *nvkat);

	// node transition matrix in sparse column format
	long *ig_t_nodes; 
	long *jg_t_nodes;
	double *gg_t_nodes;
	//edge-transition matrix in sparse column format
	long *ig_t_edges; 
	long *jg_t_edges;
	double *gg_t_edges;
	void LoadTmatrices();
	//------------------------------------------------------------------------------------------------------------------------------


	int Enumerate_Edges_In_Nonconforming_Mesh(); // numbering of edges in a nodal mesh (edges, ed generation)
	int Build_Sigma_Stucture(); // construction of ig_s, jg_s, s_val
	int Build_T_Matrix();       // generation of matrix T (final)

	// builds a chain of edge numbers needed to calculate 
	// subsequent non-zero components of column j 
	Btree* Build_Sequence(long e, double value);

	// generate a vector of weights in all edges (both terminal and non-terminal)
	int CalcValuesAll(double *v3_c, double *v3_all);
	int CalcValuesAll(double *v3); // appends to the end

	int UnloadVMesh();
};
//--------------------------------------------------------------------------------------------------