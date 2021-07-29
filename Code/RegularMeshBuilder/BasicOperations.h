#ifndef __BasicOperations__H__
#define __BasicOperations__H__
#include <set>
#include <map>

extern FILE *log_file;
using namespace std;
#define _PI_ 3.14159265358979
	int open_file_r(char *file_name, FILE **file_stream);

	int open_file_w(char *file_name, FILE **file_stream);

	int open_file_rb(char *file_name, FILE **file_stream);

	int open_file_wb(char *file_name, FILE **file_stream);

	int find_match(set<int> &a1, set<int> &a2);

	int find_match(set<int> &a1, set<int> &a2, set<int> &a3, set<int> &a4);

	int find_match(set<int> &a1, set<int> &a2, int except);

	int find_match(set<int> &a1, set<int> &a2, set<int> &a3, set<int> &a4, int except);

	template <typename T1>
	int find_match(map<int, T1> &a1, map<int, T1> &a2)
	{
		map<int, T1>::iterator it1, it2;

		for (it1 = a1.begin(); it1 != a1.end(); it1++)
		{
			if (a2.find(it1->first) != a2.end())
			{
				return it1->first;
				break;
			}
		}

		return -1;
	}

	template <typename T1>
	int find_match(map<int, T1> &a1, map<int, T1> &a2, map<int, T1> &a3, map<int, T1> &a4)
	{
		map<int, T1>::iterator it1;

		for (it1 = a1.begin(); it1 != a1.end(); it1++)
			if (a2.find(it1->first) != a2.end())
				if (a3.find(it1->first) != a3.end())
					if (a4.find(it1->first) != a4.end())
						return it1->first;
		return -1;
	}

	void M_m_V(double *M, double *V, double *R, int size);

	double V_m_V(double *V1, double *V2, int size);

	double V_m_V(double *V1, double *V2);

	void V_m_V(double *V1, double *V2, double *V3);

	double V_a_V(double *V1, double *V2, double *R);

	double V_a_V(double *V1, double *V2, double *R, int size);

	double V_s_V(double *V1, double *V2, double *R);

	double V_s_V(double *V1, double *V2, double *R, int size);

	double V_norm(double *V, int size);

	double V_norm(double *V);

	void normalize_V(double *V);

	void normalize_V(double *V, int size);

	void calc_orthogonal_vector(double *V1, double *V2, double *V3);

	void calc_normal(double *p1, double *p2, double *p3, double *n);

	void quadrangle_center(double *p1, double *p2, double *p3, double *p4, double *c);

	void quadrangle_center(double *qp, double *c);

	void element_center(int *nn, double *xyz, double *c);

	bool ray_intersects_polygon(double *p1, double *p2, double *p3, double *p4, double *r1, double *r2, double *ip);

	bool ray_intersects_triangle(double *p1, double *p2, double *p3, double *r1, double *r2, double *ip);

	double calc_distance(double *pnt1, double *pnt2);

	void determine_barycentric_coordinates(double *p1, double *p2, double *p3, double *ip, double *lc);

	double triangle_mes(double a, double b, double c);

	bool point_in_triangle(double *p1, double *p2, double *p3, double *ip);

	void calc_point_in_polygon(double &u, double &v, double *p1, double *p2, double *p3, double *p4, double *p);

	void calc_point_in_hexahedron(double &u, double &v, double &w, double *hp, double *p);

	template <typename T>
	void find_min_max(T *arr, int size, T &min_v, T &max_v)
	{
		int i;

		if (size < 1)
		{
			min_v = max_v = 0.0;
			return;
		}

		min_v = max_v = arr[0];
		for (i = 1; i < size; i++)
		{
			if (min_v > arr[i]) min_v = arr[i];
			if (max_v < arr[i]) max_v = arr[i];
		}
	}

	template <typename T>
	void find_min_max(T *arr, int size, T &min_v, T &max_v, int step)
	{
		int i;

		if (size < 1)
		{
			min_v = max_v = 0.0;
			return;
		}

		min_v = max_v = arr[0];
		for (i = step; i < size; i += step)
		{
			if (min_v > arr[i]) min_v = arr[i];
			if (max_v < arr[i]) max_v = arr[i];
		}
	}

	bool calc_local_coordinates_in_hexahedron(double *R, double *HexPnt, double *lc);

	int open_log(char *file_name);

	void write_to_log(char *str);
	void write_to_log(const char *str);
#endif