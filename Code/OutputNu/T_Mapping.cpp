#include "stdafx.h"
#include "T_Mapping.h"
#include "Portret.h"
extern ofstream logfile;
//----------------------------------------------------------------------------------------- 
int T_Mapping_Vec::UnloadVMesh()
{
	int i,j;
	int tmplong;
	double tmpdouble;
	FILE *fp;

	const int size_i=sizeof(int);
	const int size_d=sizeof(double);

	if(!(fp=fopen("nodesforedges.dat","wb"))){
		printf("Error open file nodesforedges.dat");
		return 1;
		}
	for(i=0;i<n;i++){
		tmplong=edges[i][0]+1;
		fwrite(&tmplong,size_i,1,fp);
		tmplong=edges[i][1]+1;
		fwrite(&tmplong,size_i,1,fp);
		}
	fclose(fp);

	if(!(fp=fopen("edges.dat","wb"))){
		printf("Error open file edges.dat");
		return 1;
		}
	for(i=0;i<kpar;i++){
		for(j=0;j<12;j++){
			tmplong=ed[i][j]+1;
			fwrite(&tmplong,size_i,1,fp);
			}
		tmplong=0;
		for(j=0;j<12;j++)fwrite(&tmplong,size_i,1,fp);
		tmplong=1;
		fwrite(&tmplong,size_i,1,fp);
		}
	fclose(fp);

	if(!(fp=fopen("tsize3d_.dat","w"))){
		printf("Error open file tsize3d_.dat");
		return 1;
		}
	fprintf(fp,"%d\n",n_dc);
	fprintf(fp,"%d",n);
	fclose(fp);

	if(!(fp=fopen("ig3d_.dat","wb"))){
		printf("Error open file ig3dv.dat");
		return 1;
		}
	j=n+1;
	for(i=n_c;i<j;i++){
		tmplong=ig_t[i]+1;
		fwrite(&tmplong,size_i,1,fp);
		}
	fclose(fp);

	if(!(fp=fopen("jg3d_.dat","wb"))){
		printf("Error open file jg3d_.dat");
		return 1;
		}
	j=ig_t[n];
	for(i=0;i<j;i++){
		tmplong=jg_t[i]+1;
		fwrite(&tmplong,size_i,1,fp);
		}
	fclose(fp);

	if(!(fp=fopen("gg3d_.dat","wb"))){
		printf("Error open file gg3d_.dat");
		return 1;
		}
	for(i=0;i<j;i++){
		tmpdouble=gg_t[i];
		fwrite(&tmpdouble,size_d,1,fp);
		}
	fclose(fp);
	return 0;
}
//----------------------------------------------------------------------------------------- 
T_Mapping_Vec::T_Mapping_Vec(long (*nver)[14], double (*xyz)[3], long kuzlov, long kpar)
{
	this->nver = nver;
	this->xyz = xyz;
	this->kuzlov = kuzlov;
	this->kpar = kpar;

	ed = NULL;
	edges = NULL;
	gg_t = NULL;
	ig_t = NULL;
	jg_t = NULL;
	ig_s = NULL;
	jg_s = NULL;
	s_val = NULL;
}
//-----------------------------------------------------------------------------------------
T_Mapping_Vec::~T_Mapping_Vec()
{
	if(ed) {delete [] ed; ed = NULL;}
	if(gg_t) {delete [] gg_t; gg_t = NULL;}
	if(ig_t) {delete [] ig_t; ig_t = NULL;}
	if(jg_t) {delete [] jg_t; jg_t = NULL;}
	if(ig_s) {delete [] ig_s; ig_s = NULL;}
	if(jg_s) {delete [] jg_s; jg_s = NULL;}
	if(s_val) {delete [] s_val; s_val = NULL;}
	if(edges) {delete [] edges; edges = NULL;}
}
//-------------------------------------------------------------------------
// нумерация рёбер в узловой сетке (генерация edges, ed)
//-------------------------------------------------------------------------
int T_Mapping_Vec::Enumerate_Edges_In_Nonconforming_Mesh()
{
	long type;    // тип элемента (всего 31 тип)
	list *s=NULL;      // структура для обычных рёбер
	list *s_term=NULL; // структура для терминальных рёбер
	list *l=NULL;      // для того, чтобы бегать по структуре
	Portret P;    // функции для работы со структурой хранятся в классе Portret
	long *ig_old_reg=NULL, *jg_old_reg=NULL; // с повторениями
	long *ig_new_reg=NULL, *jg_new_reg=NULL; // без повторений
	long *ig_ter=NULL, *jg_ter=NULL; // ig, jg, получаемые из структуры для терминальных рёбер
	long *orientation=NULL; // глобальная ориентация рёбер 
	long flag;
	long i, j, k, m, tmp1, tmp2, v0, v1, t;

	long reg[12][2]={ 	// какие нетерминальные рёбра есть на элементе
	1,2, 3,4, 5,6, 7,8, 1,3, 5,7, 2,4, 6,8, 1,5, 2,6, 3,7, 4,8 };

	// какие терминальные рёбра есть на элементе
	// (для шестигранника связи между локальными рёбрами и вершинами те же, что и для 
	// параллелепипеда, только тип элемента увеличен на 31, поэтому при обработке
	// шестигранника мы просто отнимаем 31 от типа эл-та).
	long ter[31][12][2]={
//	1      2      3      4      5      6      7      8      9      10     11     12
	0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 1
	1,9,   9,3,   5,10,  10,7,  9,10,  0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 2
	1,9,   9,2,   5,10,  10,6,  9,10,  0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 3
	2,9,   9,4,   6,10,  10,8,  9,10,  0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 4
	3,9,   9,4,   7,10,  10,8,  9,10,  0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 5
	1,9,   9,5,   3,10,  10,7,  9,10,  0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 6
	2,9,   9,6,   4,10,  10,8,  9,10,  0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 7
	1,9,   9,5,   2,10,  10,6,  9,10,  0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 8
	3,9,   9,7,   4,10,  10,8,  9,10,  0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 9
	1,11,  11,3,  9,13,  13,10, 5,12,  12,7,  1,9,   9,5,   11,13, 13,12, 3,10,  10,7, // type 10
	2,11,  11,4,  9,13,  13,10, 6,12,  12,8,  2,9,   9,6,   11,13, 13,12, 4,10,  10,8, // type 11
	1,11,  11,2,  9,13,  13,10, 5,12,  12,6,  1,9,   9,5,   11,13, 13,12, 2,10,  10,6, // type 12
	3,11,  11,4,  9,13,  13,10, 7,12,  12,8,  3,9,   9,7,   11,13, 13,12, 4,10,  10,8, // type 13
	1,9,   9,5,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 14
	2,9,   9,6,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 15
	4,9,   9,8,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 16
	3,9,   9,7,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 17
	2,9,   9,4,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 18
	1,9,   9,3,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 19
	3,9,   9,4,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 20
	1,9,   9,2,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 21
	1,9,   9,3,   2,10,  10,4,  9,10,  0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 22
	1,9,   9,2,   3,10,  10,4,  9,10,  0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 23
	1,11,  11,2,  9,13,  13,10, 3,12,  12,4,  1,9,   9,3,   11,13, 13,12, 2,10,  10,4, // type 24
	6,9,   9,8,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 25
	5,9,   9,7,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 26
	7,9,   9,8,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 27
	5,9,   9,6,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 28
	5,9,   9,7,   6,10,  10,8,  9,10,  0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 29
	5,9,   9,6,   7,10,  10,8,  9,10,  0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 30
	5,11,  11,6,  9,13,  13,10, 7,12,  12,8,  5,9,   9,7,   11,13, 13,12, 6,10,  10,8  // type 31  
	};

	// в Си нумерация с нуля -> вычитаем 1
	for(i=0; i<12; i++)
	for(j=0; j<2;  j++) 
		reg[i][j]--;

	for(i=0; i<31; i++)
	for(j=0; j<12; j++) 
	for(k=0; k<2;  k++)
		ter[i][j][k]--;

	s = new list[kuzlov];      // выделяем память под структуру для регулярных рёбер
	if(s == 0)
		Memory_allocation_error("s", "T_Mapping_Vec::Enumerate_Edges_In_Nonconforming_Mesh");

	s_term = new list[kuzlov]; // выделяем память под структуру для терминальных рёбер
	if(s_term == 0)
		Memory_allocation_error("s_term", "T_Mapping_Vec::Enumerate_Edges_In_Nonconforming_Mesh");

	for(i=0; i<kuzlov; i++) // обнуляем
	{
		s[i].next = NULL;
		s_term[i].next = NULL;
	}

	// цикл по эл-там (обычные рёбра заносим в одну структуру, терминальные - в другую)
	for(i=0; i<kpar; i++)
	{
		type = nver[i][13];
		if(type>30) // если шестигранник
			type -= 31;

		for(j=0; j<12; j++)
		{
			// обычные рёбра
			v0 = reg[j][0];
			v1 = reg[j][1];

			tmp1 = nver[i][v0];
			tmp2 = nver[i][v1];

			// заносим в структуру для регулярных рёбер
			if(tmp1 < tmp2) 
			{	P.Add_In_Ordered_List(&s[tmp2], tmp1);	}
			else
			{	P.Add_In_Ordered_List(&s[tmp1], tmp2);	}

			// терминальные рёбра
			v0 = ter[type][j][0];
			v1 = ter[type][j][1];

			if(v0!=-1 && v1!=-1)
			{
				tmp1 = nver[i][v0];
				tmp2 = nver[i][v1];

				// заносим в структуру для терминальных рёбер
				if(tmp1 < tmp2)
				{	P.Add_In_Ordered_List(&s_term[tmp2], tmp1);	}
				else
				{	P.Add_In_Ordered_List(&s_term[tmp1], tmp2);	}
			}
		}
	}

	// подсчитываем размер структур
	// здесь значения n_c, n_dc, n не окончательные - после удаления повторений они изменятся
	n_c  = 0; // размер s
	n_dc = 0; // размер s_term
	for(i=0; i<kuzlov; i++)
	{
		n_c += P.Size_Of_List(s[i].next);
		n_dc += P.Size_Of_List(s_term[i].next);
	}
	n = n_c + n_dc; 

	// выделяем память 
	ig_old_reg = new long[kuzlov+1];
	if(ig_old_reg == 0)
		Memory_allocation_error("ig_old_reg", "T_Mapping_Vec::Enumerate_Edges_In_Nonconforming_Mesh");

	jg_old_reg = new long[n_c];
	if(jg_old_reg == 0)
		Memory_allocation_error("jg_old_reg", "T_Mapping_Vec::Enumerate_Edges_In_Nonconforming_Mesh");

	ig_ter = new long[kuzlov+1];
	if(ig_ter == 0)
		Memory_allocation_error("ig_ter", "T_Mapping_Vec::Enumerate_Edges_In_Nonconforming_Mesh");

	jg_ter = new long[n_dc];
	if(jg_ter == 0)
		Memory_allocation_error("jg_ter", "T_Mapping_Vec::Enumerate_Edges_In_Nonconforming_Mesh");

	// заполняем ig_old, jg_old
	// обычные рёбра
	i = 0;
	ig_old_reg[0] = 0;
	for(k=0;k<kuzlov;k++)
	{
		l = s[k].next;
		while(l != NULL)
		{
			jg_old_reg[i] = l->number;
			i++;
			l = l->next;
		}
		P.Clear_List(s[k].next);
		ig_old_reg[k+1] = i;
	}

	// терминальные рёбра
	i = 0;
	ig_ter[0] = 0;
	for(k=0;k<kuzlov;k++)
	{
		l = s_term[k].next;
		while(l != NULL)
		{
			jg_ter[i] = l->number;
			i++;
			l = l->next;
		}
		P.Clear_List(s_term[k].next);
		ig_ter[k+1] = i;
	}

	// освобождаем память
	if(s) {delete [] s; s = NULL;}
	if(s_term) {delete [] s_term; s_term = NULL;}

	// удаляем повторения (терминальные рёбра занесены 2 раза - дублируются в структуре для обычных рёбер)
	for(i=0; i<kuzlov; i++)
	{
		for(j=ig_ter[i]; j<=ig_ter[i+1]-1; j++)
		{
			k = jg_ter[j];
			// теперь ищем ребро (i,k) в структуре для обычных рёбер
			for(m=ig_old_reg[i]; m<=ig_old_reg[i+1]-1; m++)
				if(jg_old_reg[m] == k)
				{
					jg_old_reg[m] = -2;
					break;
				}
		}
	}

	// переписываем ребра в массив edges[][2] (без повторений)...
	edges = new long[n_c][2]; 
	if(edges == 0)
		Memory_allocation_error("edges", "T_Mapping_Vec::Enumerate_Edges_In_Nonconforming_Mesh");

	k = 0;
    for(i=0; i<kuzlov; i++) // обычные рёбра
	{
		for(j=ig_old_reg[i]; j<=ig_old_reg[i+1]-1; j++)
			if(jg_old_reg[j] >= 0)
			{
				edges[k][1] = i;
				edges[k][0] = jg_old_reg[j];
				k++;
			}
	}

    for(i=0; i<kuzlov; i++) // терминальные рёбра
	{
		for(j=ig_ter[i]; j<=ig_ter[i+1]-1; j++)
		{
			edges[k][1] = i;
			edges[k][0] = jg_ter[j];
			k++;
		}
	}

	// истинное число обычных и терминальных рёбер (без повторений)
	n = n_c;
	n_c = n_c - n_dc;
	
	// новые ig_new_reg[], jg_new_reg[]
	ig_new_reg = new long[kuzlov+1];
	if(ig_new_reg == 0)
		Memory_allocation_error("ig_new_reg", "T_Mapping_Vec::Enumerate_Edges_In_Nonconforming_Mesh");

	jg_new_reg = new long[n_c];
	if(jg_new_reg == 0)
		Memory_allocation_error("jg_new_reg", "T_Mapping_Vec::Enumerate_Edges_In_Nonconforming_Mesh");

	ig_new_reg[0] = 0;
	k = 0;
	for(i=0; i<kuzlov; i++)
	{
		for(j=ig_old_reg[i]; j<=ig_old_reg[i+1]-1; j++)
		{
			if(jg_old_reg[j] >= 0)
			{
				jg_new_reg[k] = jg_old_reg[j];
				k++;
			}
		}
		ig_new_reg[i+1] = k;
	}

	// освобождаем память
	if(ig_old_reg) {delete [] ig_old_reg; ig_old_reg=NULL;}
	if(jg_old_reg) {delete [] jg_old_reg; jg_old_reg=NULL;}

	// теперь нужно сгенерировать структуру, в которой для каждого элемента храняться номера его рёбер + 
	// номера терминальных рёбер (ed[][25]).
	ed = new long[kpar][25];
	if(ed == 0)
		Memory_allocation_error("ed", "T_Mapping_Vec::Enumerate_Edges_In_Nonconforming_Mesh");

	for(i=0; i<kpar; i++)
	{
		type = nver[i][13]; // тип		
		if(type>30)
			type -= 31;

		ed[i][24] = type; 
		for(j=0; j<12; j++) 
		{
			// обычные рёбра
			v0 = nver[i][reg[j][0]];
			v1 = nver[i][reg[j][1]];
			Sort2(&v0, &v1);
			flag = 0;
			for(k=ig_new_reg[v0]; k<=ig_new_reg[v0+1]-1; k++)
				if(jg_new_reg[k] == v1)
				{
					ed[i][j] = k;
					flag = 1;
					break;
				}
			// если не нашли ребро(v0,v1) в структуре для регулярных рёбер -
			// значит оно терминальное
			if(flag == 0) 
			{
				for(k=ig_ter[v0]; k<=ig_ter[v0+1]-1; k++)
					if(jg_ter[k] == v1)
					{
						ed[i][j] = k + n_c;
						break;
					}
			}

			// терминальные рёбра
			tmp1 = ter[type][j][0];
			tmp2 = ter[type][j][1];
			if(tmp1 >=0 && tmp2>=0)
			{
				v0 = nver[i][tmp1];
				v1 = nver[i][tmp2];
				Sort2(&v0, &v1);
				for(k=ig_ter[v0]; k<=ig_ter[v0+1]-1; k++)
					if(jg_ter[k] == v1)
					{
						ed[i][j+12] = k + n_c;
						break;
					}
			}
			else
			{
				ed[i][j+12] = -2;
			}
		}
	}

	// освобождаем память
	if(ig_ter) {delete [] ig_ter; ig_ter=NULL;}
	if(jg_ter) {delete [] jg_ter; jg_ter=NULL;}
	if(ig_new_reg) {delete [] ig_new_reg; ig_new_reg=NULL;}
	if(jg_new_reg) {delete [] jg_new_reg; jg_new_reg=NULL;}

	// ориентация (определяем глобальную ориентацию по локальной ориентации)
	orientation = new long[n];	
	if(orientation == 0)
		Memory_allocation_error("orientation", "T_Mapping_Vec::Enumerate_Edges_In_Nonconforming_Mesh");

	for(i=0; i<this->n; i++)
		orientation[i] = 1;

	for(i=0; i<this->kpar; i++)
	{
		long v[8];

		for(j=0; j<8; j++)
			v[j] = nver[i][j];

		for(j=0; j<12; j++)
			if(v[reg[j][0]] > v[reg[j][1]]) 
			{// глобальная ориентация не соответствует локальной
				orientation[ed[i][j]] = -1;
			}
	}
	// если глобальная ориентация не соответствует локальной, 
	// то меняем местами начало и конец ребра в edges[][2]
	for(i=0; i<n; i++)
	{
		if (orientation[i]==-1)
		{
			t = edges[i][0];
			edges[i][0] = edges[i][1];
			edges[i][1] = t;
		}
	}

	// освобождаем память
	if(orientation) {delete [] orientation; orientation=NULL;}

	return 0;
}
//-----------------------------------------------------------------------------------
//- Построение вспомогательной структуры для Т-матрицы создание ig_s, jg_s, s_val
//-----------------------------------------------------------------------------------
int T_Mapping_Vec::Build_Sigma_Stucture()
{
	// структура sigma хранится в виде массивов ig_s и jg_s
	list *s=NULL;       // структура
	long size_s;   
	list *l=NULL;       // для того, чтобы бегать по структуре
	Portret P;     // функции для работы со структурой хранятся в классе Portret
	long i, j, k, m, e, e1, e2;
	long type;
	In_Out R;
	double x0_reg[3], x1_reg[3], axis_reg[3];
	double x0_ter[3], x1_ter[3], axis_ter[3];
	double x0_reg2[3];
	double a, b, len;

    // локальные номера рёбер КЭ (содержащего терминальное ребро),которым соответствуют локальные базисные функции
	// принимающие ненулевые значения в рассматриваемом терминальном ребре
		long n_sig_ed[31][12][2] = {
//	1      2      3      4      5      6      7      8      9      10     11     12
	0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 1
	5,0,   5,0,   6,0,   6,0,   9,11,  0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 2
	1,0,   1,0,   3,0,   3,0,   9,10,  0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 3
	7,0,   7,0,   8,0,   8,0,   10,12, 0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 4
	2,0,   2,0,   4,0,   4,0,   11,12, 0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 5
	9,0,   9,0,   11,0,  11,0,  5,6,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 6
	10,0,  10,0,  12,0,  12,0,  7,8,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 7
	9,0,   9,0,   10,0,  10,0,  1,3,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 8
	11,0,  11,0,  12,0,  12,0,  2,4,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 9
	5,0,   5,0,   13,17, 14,18, 6,0,   6,0,   9,0,   9,0,   19,23, 20,24, 11,0,  11,0, // type 10
	7,0,   7,0,   13,17, 14,18, 8,0,   8,0,   10,0,  10,0,  19,23, 20,24, 12,0,  12,0, // type 11
	1,0,   1,0,   13,17, 14,18, 3,0,   3,0,   9,0,   9,0,   19,23, 20,24, 10,0,  10,0, // type 12
	2,0,   2,0,   13,17, 14,18, 4,0,   4,0,   11,0,  11,0,  19,23, 20,24, 12,0,  12,0, // type 13
	9,0,   9,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 14
	10,0,  10,0,  0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 15
	12,0,  12,0,  0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 16
	11,0,  11,0,  0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 17
	7,0,   7,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 18
	5,0,   5,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 19
	2,0,   2,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 20
	1,0,   1,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 21
	5,0,   5,0,   7,0,   7,0,   1,2,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 22
	1,0,   1,0,   2,0,   2,0,   5,7,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 23
	1,0,   1,0,   13,17, 14,18, 2,0,   2,0,   5,0,   5,0,   19,23, 20,24, 7,0,   7,0,  // type 24
	8,0,   8,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 25
	6,0,   6,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 26
	4,0,   4,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 27
	3,0,   3,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 28
	6,0,   6,0,   8,0,   8,0,   3,4,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 29
	3,0,   3,0,   4,0,   4,0,   6,8,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  // type 30
	3,0,   3,0,   13,17, 14,18, 4,0,   4,0,   6,0,   6,0,   19,23, 20,24, 8,0,   8,0
	};

	 // в Си нумерация с нуля -> вычитаем 1
	for(i=0;i<31;i++)
	for(j=0;j<12;j++)
	for(k=0;k<2;k++)
		n_sig_ed[i][j][k]--;

	s = new list[n]; // выделяем память под структуру
	if(s == 0)
		Memory_allocation_error("s", "T_Mapping_Vec::Build_Sigma_Stucture");

	for(i=0;i<n;i++)
		s[i].next = NULL;

	for(i=0; i<kpar; i++)
	{
		type = ed[i][24];
		for(j=0; j<12; j++) // максимум 12 терминальных рёбер на элементе
			for(k=0; k<2; k++)
			{
				e = n_sig_ed[type][j][k];
				if(e != -1)
				{
					P.Add_In_Ordered_List(&s[ed[i][j+12]], ed[i][e]);
				}
			}//k
	}

	// подсчитываем размер структуры
	size_s = 0;
	for(i=0; i<n; i++)
		size_s += P.Size_Of_List(s[i].next);

	// выделяем память
	ig_s = new long[n+1];
	if(ig_s == 0)
		Memory_allocation_error("ig_s", "T_Mapping_Vec::Build_Sigma_Stucture");

	jg_s = new long[size_s];
	if(jg_s == 0)
		Memory_allocation_error("jg_s", "T_Mapping_Vec::Build_Sigma_Stucture");

	// заполняем ig_s, jg_s
	i = 0;
	ig_s[0] = 0;
	for(k=0; k<n; k++)
	{
		l = s[k].next;
		while(l != NULL)
		{
			this->jg_s[i] = l->number;
			i++;
			l = l->next;
		}
		P.Clear_List(s[k].next);
		this->ig_s[k+1] = i;
	}

	if(s) {delete [] s; s=NULL;}

	// заполняем s_val
	s_val = new double[size_s];
	if(s_val == 0)
		Memory_allocation_error("s_val", "T_Mapping_Vec::Build_Sigma_Stucture");

	for(i=0; i<kpar; i++)
	{
		type = ed[i][24];
		for(j=0; j<12; j++) // максимум 12 терминальных рёбер на элементе
		{
			e1 = n_sig_ed[type][j][0]; // локальные 
			e2 = n_sig_ed[type][j][1];

			// на ребре лежит терминальное ребро
			if(e1 !=-1 && e2 == -1) 
			{                      
				e = ed[i][j+12]; // терминальное ребро (глобальный номер)
				e1 = ed[i][e1];  // ребро, к-рое содержит терминальное ребро

				for(k=0;k<3;k++)
				{
					x0_reg[k] = xyz[edges[e1][0]][k];
					x1_reg[k] = xyz[edges[e1][1]][k];

					x0_ter[k] = xyz[edges[e][0]][k];
					x1_ter[k] = xyz[edges[e][1]][k];

					axis_reg[k] = x1_reg[k] - x0_reg[k];
					axis_ter[k] = x1_ter[k] - x0_ter[k];
				}

				for(m=ig_s[e]; m<=ig_s[e+1]-1; m++)
				{
					if(jg_s[m]==e1)
					{
						s_val[m] = Norm_Euclid(axis_ter, 3)/Norm_Euclid(axis_reg, 3);
						break;
					}
				}
			}

			// терминальное ребро лежит в грани)
			else if(e1 != -1 && e2 != -1) // если сразу 2 функции принимают ненулевое значение на терминальном ребре
			{
				// теперь глобальные
				e = ed[i][j+12]; // терминальное ребро
				e1 = ed[i][e1];  // два нетерминальных(?) ребра, участвующих в построении функции на терминальном ребре
				e2 = ed[i][e2];

				for(k=0; k<3; k++)
				{
					x0_ter[k]  = xyz[edges[e][0]][k];
					x0_reg[k]  = xyz[edges[e1][0]][k];
					x0_reg2[k] = xyz[edges[e2][0]][k];
				}

				a = Interval(x0_ter, x0_reg);
				b = Interval(x0_ter, x0_reg2);
				len = a + b;

				for(m=ig_s[e]; m<=ig_s[e+1]-1; m++)
				{
					if(jg_s[m]==e1)
					{
						s_val[m] =	b/len;
					}
					else if(jg_s[m]==e2)
					{
						s_val[m] = a/len;
					}
				}
			}	
		}//j
	}

	return 0;
}
//-------------------------------------------------------------------------------------------------
//--- генерация матрицы T (окончательная)
//--------------------------------------------------------------------------------------------------
int T_Mapping_Vec::Build_T_Matrix()
{
	long i, j, k;
	std::vector<long_double> *s_t=NULL; // структура (для матрицы T; массивы ig_t, jg_t, gg_t создаются из неё)
	long_double tmp;
	Btree *d=NULL, *dd=NULL;
	std::stack<Btree*> s;
	long size_s_t; 

	s_t = new std::vector<long_double>[n];
	if(s_t == 0)
		Memory_allocation_error("s_t", "T_Mapping_Vec::Build_T_Matrix");

	for(j=this->n_c; j<this->n; j++) // Выбираем очередной номер столбца j (начиная с j=n_c+1). 
	{
		// Для него определяем номера рёбер k и соответствующие им значения $T_{kj}$  или $\tilda T_{kj}$
		// по сформированной ранее структуре данных Sigma.

		for(i=this->ig_s[j]; i<=this->ig_s[j+1]-1; i++)
		{
			k = jg_s[i];
			if(k < this->n_c) // просто добавляем
			{ 
				tmp.i = k;
				tmp.d = s_val[i];
				s_t[j].push_back(tmp);
			}
			else // k>=n_c
			{ // начинается построение цепочки
				dd = this->Build_Sequence(k, s_val[i]); // построили бинарное дерево

				// обходим дерево и заносим соответствующие значения
				d = dd;
				d->visit(s_t, d, j);
				delete dd; dd=NULL;
			}//else  k>=n_c
		}// i
	} // j

	// освобождаем память для вспомогательной структуры 
	if(s_val) {delete [] s_val;	s_val = NULL;}
	if(ig_s)  {delete [] ig_s;	ig_s = NULL;}
	if(jg_s)  {delete [] jg_s;	jg_s = NULL;}

	// подсчитываем размер структуры
	size_s_t = 0;
	for(i=0; i<n; i++)
		size_s_t += (long)s_t[i].size();

	// выделяем память под T-матрицу
	ig_t = new long[n+1];
	if(ig_t == 0)
		Memory_allocation_error("ig_t", "T_Mapping_Vec::Build_T_Matrix");

	jg_t = new long[size_s_t];
	if(jg_t == 0)
		Memory_allocation_error("jg_t", "T_Mapping_Vec::Build_T_Matrix");

	gg_t = new double[size_s_t];
	if(gg_t == 0)
		Memory_allocation_error("gg_t", "T_Mapping_Vec::Build_T_Matrix");

	// заполняем ig_t, jg_t, gg_t
	i = 0;
	this->ig_t[0] = 0;
	for(k=0; k<n; k++)
	{
		for(j=0; j<(long)s_t[k].size(); j++)
		{
			this->jg_t[i] = s_t[k][j].i;
			this->gg_t[i] = s_t[k][j].d;
			i++;
		}
		this->ig_t[k+1] = i;
	}

	// удаляем s_t
	for(i=0;i<n;i++)
	{
		long size;
		size = (long)s_t[i].size();
		for(j=0;j<size;j++)
			s_t[i].pop_back();
	}

	if(s_t) {delete [] s_t; s_t = NULL;}

	FILE *t_size=NULL;
	t_size = fopen("t_size.dat", "w");
	if(t_size==0)
	{
		Cannot_open_file("t_size.dat", "T_Mapping_Vec::Build_T_Matrix");
	}
	else
	{
		fprintf(t_size, "%ld\t terminal edges\n", n_dc);
		fprintf(t_size, "%ld\t regular edges\n", n_c);
		fprintf(t_size, "%ld\t edges total\n", n);
		fclose(t_size);
	}

	return 0;
}
//-----------------------------------------------------------------
// строит цепочку номеров рёбер, необходимых для вычисления 
// последующих ненулевых компонент столбца j
//-----------------------------------------------------------------
Btree* T_Mapping_Vec::Build_Sequence(long e, double value)
{
	Btree *d=NULL;
	long l, r;
	long n;
	double new_val;
	
	n = ig_s[e+1] - ig_s[e]; // число функций, принимающих ненулевые значения
	// в ребре с глобальным номером e

	d = new Btree(e, value);
	if(d == 0)
		Memory_allocation_error("d", "T_Mapping_Vec::Build_Sequence");

	if(n==1)
	{
		l = jg_s[ig_s[e]];
		new_val = value * s_val[ig_s[e]];
		d->left = Build_Sequence(l, new_val);
		return d;
	}
	else if(n==2)
	{
		l = jg_s[ig_s[e]];
		new_val = value * s_val[ig_s[e]];
		d->left = Build_Sequence(l, new_val);

		r = jg_s[ig_s[e]+1];
		new_val = value * s_val[ig_s[e]+1];
		d->right = Build_Sequence(r, new_val);

		return d;
	}
	else 
	{
		return d;
	}
}
//------------------------------------------------------------------------
// сформировать вектор весов во всех рёбрах (и в терминальных, и в нетерминальных)
//------------------------------------------------------------------------
int T_Mapping_Vec::CalcValuesAll(double *v3_c, double *v3_all)
{
	long i, j;
	double s;

	for (i=0; i<n_c; i++)
		v3_all[i] = v3_c[i];

	for (i=n_c; i<n; i++)
	{
		s = 0.0;
		for (j=ig_t[i]; j<=ig_t[i+1]-1; j++)
			s += v3_c[jg_t[j]]*gg_t[j];
		v3_all[i] = s;
	}
	return 0;
}
//------------------------------------------------------------------------
// сформировать вектор весов во всех рёбрах (и в терминальных, и в нетерминальных)
// (дописывает в конец того же вектора)
//------------------------------------------------------------------------
int T_Mapping_Vec::CalcValuesAll(double *v3)
{
	long i, j;
	double s;

	for (i=n_c; i<n; i++)
	{
		s = 0.0;
		for (j=ig_t[i]; j<=ig_t[i+1]-1; j++)
			s += v3[jg_t[j]]*gg_t[j];
		v3[i] = s;
	}
	return 0;
}
//------------------------------------------------------------------
Btree::Btree()
{
	this->left = 0;
	this->right = 0;
	this->visited = false;
}
//-------------------------------------------------------------------
Btree::Btree(long elem_long, double elem_double)
{
	this->elem_double = elem_double;
	this->elem_long = elem_long;
	this->left = 0;
	this->right = 0;	
	this->visited = false;
}
//--------------------------------------------------------------------
Btree::~Btree()
{
	if(this->left != 0)
	{
		DeleteTree(this->left);
		this->left = 0;
	}
	if(this->right != 0)
	{
		DeleteTree(this->right);
		this->right = 0;
	}
}
//--------------------------------------------------------------------
int Btree::DeleteTree(Btree *t)
{
	if(t!=0)
	{
		DeleteTree(t->left);
		t->left = 0;
		DeleteTree(t->right);
		t->right = 0;
		delete t; t=NULL;
		return 0;
	}
	return 0;
}
//---------------------------------------------------------------------
Btree *Btree::visit(std::vector<long_double> *s_t, Btree *t, long j)
{
	long_double tmp;

	if(t!=0)
	{
		visit(s_t,t->left,j);
		visit(s_t,t->right,j);

		if(t->left == 0 && t->right ==0)
		{
			tmp.i = t->elem_long;
			tmp.d = t->elem_double;
			s_t[j].push_back(tmp);
		}		
	}
	return 0;	
}
//---------------------------------------------------------------------
int Btree::Add_Left(long elem_long, double elem_double)
{
	if(this->left!=0)
	{
		cerr << "Error in function Btree::Add_Left: left subtree != NULL\n";
		logfile << "Error in function Btree::Add_Left: left subtree != NULL\n";
		throw logic_error("Error in function Btree::Add_Left: left subtree != NULL\n");
		return 1;
	}

	this->left = new Btree(elem_long, elem_double);
	if(this->left == 0)
		Memory_allocation_error("this->left", "Btree::Add_Left");

	return 0;
}
//----------------------------------------------------------------------
int Btree::Add_Right(long elem_long, double elem_double)
{
	if(this->right!=0)
	{
		cerr << "Error in function Btree::Add_Right: right subtree != NULL\n";
		logfile << "Error in function Btree::Add_Right: right subtree != NULL\n";
        throw logic_error("Error in function Btree::Add_Right: right subtree != NULL\n");
		return 1;
	}

	this->right = new Btree(elem_long, elem_double);
	if(this->right == 0)
		Memory_allocation_error("this->right", "Btree::Add_Right");

	return 0;
}
//-----------------------------------------------------------------------