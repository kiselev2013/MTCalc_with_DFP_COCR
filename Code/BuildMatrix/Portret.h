#pragma once
//----------------------------
struct list // список целых чисел
{
	long number;
	list *next;
};
//----------------------------
class Portret
{
public:

	long *ig;
	long *jg;
	list *s; // структура
	bool is_mem_ig_jg_allocated;
	long n; // размерность глобальной матрицы
	long n_of_elements;
	long size_jg;
	long (*nvtr)[8];
	long (*nvtr4)[4];
	long (*nver)[14];

	long n_c; // число нетерминальных узлов
	long *ig_t; // матрица T-преобразования
	long *jg_t;
	
	Portret(long (*nver)[14], long n_of_elements, long n_of_nodes, long n_c, long *ig_t, long *jg_t);
	Portret(long (*nvtr)[4],  long n_of_elements, long n_of_nodes, long n_c, long *ig_t, long *jg_t); // for 2d
	Portret(long (*nvtr)[8],  long n_of_elements, long n_of_nodes, long n_c, long *ig_t, long *jg_t);
	Portret(long (*nvtr4)[4], long n_of_elements, long n_of_nodes);
	
	Portret();
	~Portret();

	void Gen_T_Portrait(); // для 3D
	void Gen_T_Portrait2(); // для 3D
	void Gen_Portret_2d_Linear(); // для сплайна
	void Gen_Portret_2d_rect_T(); // для 2D

	void Add_In_Ordered_List(list *s, long x); // вставка эл-та в упорядоченный список
	void Clear_List(list *s);
	long Size_Of_List(list *s);
};
