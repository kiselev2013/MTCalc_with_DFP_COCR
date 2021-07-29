#include "BasicOperations.h"

FILE *log_file;

	int open_file_r(char *file_name, FILE **file_stream)
	{
		char buf[2048];

		if (!((*file_stream) = fopen(file_name, "r")))
		{
			sprintf(buf, "Error : Could not read file '%s'\n", file_name);
			write_to_log(buf);
			return 1;
		}

		return 0;
	}

	int open_file_w(char *file_name, FILE **file_stream)
	{
		char buf[2048];

		if (!((*file_stream) = fopen(file_name, "w")))
		{
			sprintf(buf, "Error : Could not write file '%s'\n", file_name);
			write_to_log(buf);
			return 1;
		}

		return 0;
	}

	int open_file_rb(char *file_name, FILE **file_stream)
	{
		char buf[2048];

		if (!((*file_stream) = fopen(file_name, "rb")))
		{
			sprintf(buf, "Error : Could not read file '%s'\n", file_name);
			write_to_log(buf);
			return 1;
		}

		return 0;
	}

	int open_file_wb(char *file_name, FILE **file_stream)
	{
		char buf[2048];

		if (!((*file_stream) = fopen(file_name, "wb")))
		{
			sprintf(buf, "Error : Could not read write '%s'\n", file_name);
			write_to_log(buf);
			return 1;
		}

		return 0;
	}

	int find_match(set<int> &a1, set<int> &a2)
	{
		set<int>::iterator it;

		for (it = a1.begin(); it != a1.end(); ++it)
			if (a2.find(*it) != a2.end())
				return *it;

		return -1;
	}

	int find_match(set<int> &a1, set<int> &a2, set<int> &a3, set<int> &a4)
	{
		set<int>::iterator it;

		for (it = a1.begin(); it != a1.end(); ++it)
			if (a2.find(*it) != a2.end())
				if (a3.find(*it) != a3.end())
					if (a4.find(*it) != a4.end())
						return *it;

		return -1;
	}

	int find_match(set<int> &a1, set<int> &a2, int except)
	{
		set<int>::iterator it;

		for (it = a1.begin(); it != a1.end(); ++it)
			if (*it != except)
				if (a2.find(*it) != a2.end())
					return *it;

		return -1;
	}

	int find_match(set<int> &a1, set<int> &a2, set<int> &a3, set<int> &a4, int except)
	{
		set<int>::iterator it;

		for (it = a1.begin(); it != a1.end(); ++it)
			if (*it != except)
				if (a2.find(*it) != a2.end())
					if (a3.find(*it) != a3.end())
						if (a4.find(*it) != a4.end())
							return *it;

		return -1;
	}

	void M_m_V(double *M, double *V, double *R, int size)
	{
		int i, j;
		double VV[3];

		VV[0] = V[0];
		VV[1] = V[1];
		VV[2] = V[2];

		for (i = 0; i < size; i++)
		{
			R[i] = M[i*size] * VV[0];
			for (j = 1; j < size; j++)
				R[i] += M[i*size + j] * VV[j];
		}
	}

	double V_m_V(double *V1, double *V2, int size)
	{
		int i;
		double res;

		res = 0.0;
		for (i = 0; i < size; i++)
			res += V1[i] * V2[i];

		return res;
	}

	double V_m_V(double *V1, double *V2)
	{
		int i;
		double res;

		res = 0.0;
		for (i = 0; i < 3; i++)
			res += V1[i] * V2[i];

		return res;
	}

	void V_m_V(double *V1, double *V2, double *V3)
	{
		V3[0] = V1[1] * V2[2] - V1[2] * V2[1];
		V3[1] = V1[2] * V2[0] - V1[0] * V2[2];
		V3[2] = V1[0] * V2[1] - V1[1] * V2[0];
	}

	double V_a_V(double *V1, double *V2, double *R)
	{
		int i;
		double res;

		for (i = 0; i < 3; i++)
			R[i] = V1[i] + V2[i];

		return res;
	}

	double V_a_V(double *V1, double *V2, double *R, int size)
	{
		int i;
		double res;

		for (i = 0; i < size; i++)
			R[i] = V1[i] + V2[i];

		return res;
	}

	double V_s_V(double *V1, double *V2, double *R)
	{
		int i;
		double res;

		for (i = 0; i < 3; i++)
			R[i] = V1[i] - V2[i];

		return res;
	}

	double V_s_V(double *V1, double *V2, double *R, int size)
	{
		int i;
		double res;

		for (i = 0; i < size; i++)
			R[i] = V1[i] - V2[i];

		return res;
	}

	double V_norm(double *V, int size)
	{
		return sqrt(V_m_V(V, V, size));
	}

	double V_norm(double *V)
	{
		return sqrt(V_m_V(V, V));
	}

	void normalize_V(double *V)
	{
		int i;
		double norm;

		norm = sqrt(V_m_V(V, V));

		for (i = 0; i < 3; i++)
			V[i] /= norm;
	}

	void normalize_V(double *V, int size)
	{
		int i;
		double norm;

		norm = sqrt(V_m_V(V, V, size));

		for (i = 0; i < size; i++)
			V[i] /= norm;
	}

	void calc_orthogonal_vector(double *V1, double *V2, double *V3)
	{
		V3[0] = -V1[1] * V2[2] + V1[2] * V2[1];
		V3[1] = -V1[2] * V2[0] + V1[0] * V2[2];
		V3[2] = -V1[0] * V2[1] + V1[1] * V2[0];
	}

	void calc_normal(double *p1, double *p2, double *p3, double *n)
	{
		double v1[9], v2[9], colinear[3];

		V_s_V(p2, p1, v1);
		V_s_V(p3, p2, v2);
		colinear[0] = fabs(V_m_V(v1, v2));

		V_s_V(p2, p1, v1 + 3);
		V_s_V(p3, p1, v2 + 3);
		colinear[1] = fabs(V_m_V(v1 + 3, v2 + 3));

		V_s_V(p3, p1, v1 + 6);
		V_s_V(p3, p2, v2 + 6);
		colinear[2] = fabs(V_m_V(v1 + 6, v2 + 6));

		if (colinear[0] < colinear[1] && colinear[0] < colinear[2])
			V_m_V(v1, v2, n);
		else
			if (colinear[1] < colinear[2])
				V_m_V(v1 + 3, v2 + 3, n);
			else
				V_m_V(v1 + 6, v2 + 6, n);

		normalize_V(n);
	}

	void quadrangle_center(double *p1, double *p2, double *p3, double *p4, double *c)
	{
		for (int i = 0; i < 3; i++)
			c[i] = (p1[i] + p2[i] + p3[i] + p4[i])*0.25;
	}

	void quadrangle_center(double *qp, double *c)
	{
		for (int i = 0; i < 3; i++)
		{
			c[i] = 0.0;
			for (int j = 0; j < 4; j++)		c[i] += qp[i + j * 3];
			c[i] *= 0.25;
		}
	}

	void element_center(int *nn, double *xyz, double *c)
	{
		int i, j, k;

		for (k=0; k<3; k++)			c[k] = xyz[nn[0]*3+k];
		for (i=1; i<8; i++)
			for (k=0; k<3; k++)			c[k] += xyz[nn[i]*3+k];

		for (k=0; k<3; k++)
			c[k] *= 0.125;
	}

	bool ray_intersects_polygon(double *p1, double *p2, double *p3, double *p4, double *r1, double *r2, double *ip)
	{
		if (ray_intersects_triangle(p1, p2, p3, r1, r2, ip) == true)	return true;
		if (ray_intersects_triangle(p1, p3, p4, r1, r2, ip) == true)	return true;

		return false;
	}

	bool ray_intersects_triangle(double *p1, double *p2, double *p3, double *r1, double *r2, double *ip)
	{

		int j;
		double d, a1, a2, a3, n[3], dir[3], ksi, ang, den, v1[3], v2[3], v3[3];

		for (j = 0; j < 3; j++)		dir[j] = r2[j] - r1[j];

		n[0] = (p2[1] - p1[1])*(p3[2] - p1[2]) - (p2[2] - p1[2])*(p3[1] - p1[1]);
		n[1] = (p2[2] - p1[2])*(p3[0] - p1[0]) - (p2[0] - p1[0])*(p3[2] - p1[2]);
		n[2] = (p2[0] - p1[0])*(p3[1] - p1[1]) - (p2[1] - p1[1])*(p3[0] - p1[0]);
		normalize_V(n);

		d = -n[0] * p1[0] - n[1] * p1[1] - n[2] * p1[2];

		den = V_m_V(n, dir);
		if (fabs(den) < 1e-10)	return false;

		ksi = -(d + n[0] * r1[0] + n[1] * r1[1] + n[2] * r1[2]) / den;

		if (ksi < 0.0 || ksi > 1.0)	return false;

		ip[0] = r1[0] + ksi * (r2[0] - r1[0]);
		ip[1] = r1[1] + ksi * (r2[1] - r1[1]);
		ip[2] = r1[2] + ksi * (r2[2] - r1[2]);

		for (j = 0; j < 3; j++)		v1[j] = p1[j] - ip[j];
		for (j = 0; j < 3; j++)		v2[j] = p2[j] - ip[j];
		for (j = 0; j < 3; j++)		v3[j] = p3[j] - ip[j];

		if (V_norm(v1) < 1e-10)
			return true;
		if (V_norm(v2) < 1e-10)
			return true;
		if (V_norm(v3) < 1e-10)
			return true;

		normalize_V(v1);
		normalize_V(v2);
		normalize_V(v3);

		a1 = V_m_V(v1, v2);
		a2 = V_m_V(v2, v3);
		a3 = V_m_V(v1, v3);

		if (a1 > 1.0) a1 = 1.0 - 1e-13;
		if (a2 > 1.0) a2 = 1.0 - 1e-13;
		if (a3 > 1.0) a3 = 1.0 - 1e-13;

		if (a1 < -1.0) a1 = -1.0 + 1e-13;
		if (a2 < -1.0) a2 = -1.0 + 1e-13;
		if (a3 < -1.0) a3 = -1.0 + 1e-13;

		ang = (acos(a1) + acos(a2) + acos(a3)) * 180.0 / _PI_;
		if (fabs(ang - 360) > 1e-4)	return false;

		return true;
	}

	double calc_distance(double *pnt1, double *pnt2)
	{
		int k;
		double d;

		for (d = 0.0, k = 0; k < 3; k++)
			d += (pnt1[k] - pnt2[k])*(pnt1[k] - pnt2[k]);
		return sqrt(d);
	}

	void determine_barycentric_coordinates(double *p1, double *p2, double *p3, double *ip, double *lc)
	{
		double d1, d2, d3, dd1, dd2, dd3, s, s1, s2, s3;

		d1 = calc_distance(p1, p2);
		d2 = calc_distance(p2, p3);
		d3 = calc_distance(p3, p1);

		dd1 = calc_distance(p1, ip);
		dd2 = calc_distance(p2, ip);
		dd3 = calc_distance(p3, ip);

		s = triangle_mes(d1, d2, d3);
		s1 = triangle_mes(d1, dd1, dd2);
		s2 = triangle_mes(d2, dd2, dd3);
		s3 = triangle_mes(d3, dd3, dd1);

		lc[0] = s2 / s;
		lc[1] = s3 / s;
		lc[2] = s1 / s;
	}

	double triangle_mes(double a, double b, double c)
	{
		double p, a1, a2, a3;
		p = 0.5*(a + b + c);
		a1 = p - a;	if (a1 < 0.0)	return 0.0;
		a2 = p - b;	if (a2 < 0.0)	return 0.0;
		a3 = p - c;	if (a3 < 0.0)	return 0.0;
		return sqrt(p*a1*a2*a3);
	}

	bool point_in_triangle(double *p1, double *p2, double *p3, double *ip)
	{
		double lc[3];

		determine_barycentric_coordinates(p1, p2, p3, ip, lc);
		if (fabs(lc[0] + lc[1] + lc[2] - 1.0) < 1e-4)
			return true;

		return false;
	}

	void calc_point_in_polygon(double &u, double &v, double *p1, double *p2, double *p3, double *p4, double *p)
	{
		if (u < 0.0)	u = 0.0;
		if (u > 1.0)	u = 1.0;
		if (v < 0.0)	v = 0.0;
		if (v > 1.0)	v = 1.0;

		for (int i = 0; i < 3; i++)
			p[i] = (1.0 - u)*(1.0 - v)*p1[i] + u*(1.0 - v)*p2[i] + u*v*p3[i] + (1.0 - u)*v*p4[i];
	}

	void calc_point_in_hexahedron(double &u, double &v, double &w, double *hp, double *p)
	{
		double q[8];
		if (u < 0.0)	u = 0.0;
		if (u > 1.0)	u = 1.0;
		if (v < 0.0)	v = 0.0;
		if (v > 1.0)	v = 1.0;
		if (w < 0.0)	w = 0.0;
		if (w > 1.0)	w = 1.0;

		q[0] = (1.0 - u)*(1.0 - v)*(1.0 - w);
		q[1] = (u)*(1.0 - v)*(1.0 - w);
		q[2] = (1.0 - u)*(v)*(1.0 - w);
		q[3] = (u)*(v)*(1.0 - w);
		q[4] = (1.0 - u)*(1.0 - v)*(w);
		q[5] = (u)*(1.0 - v)*(w);
		q[6] = (1.0 - u)*(v)*(w);
		q[7] = (u)*(v)*(w);

		for (int i = 0; i < 3; i++)
		{
			p[i] = 0.0;
			for (int k = 0; k < 8; k++)
				p[i] += q[k] * hp[k * 3 + i];
		}
	}

	void GetGlobalCoordinates(double *HexPnt, double *lc, double *gc)
	{
		double ksi, eta, phi;

		ksi = lc[0];
		eta = lc[1];
		phi = lc[2];

		gc[0] = HexPnt[0 * 3] * (1 - ksi)*(1 - eta)*(1 - phi) + HexPnt[1 * 3] * (ksi)*(1 - eta)*(1 - phi) +
			HexPnt[2 * 3] * (1 - ksi)*(eta)*(1 - phi) + HexPnt[3 * 3] * (ksi)*(eta)*(1 - phi) +
			HexPnt[4 * 3] * (1 - ksi)*(1 - eta)*(phi)+HexPnt[5 * 3] * (ksi)*(1 - eta)*(phi)+
			HexPnt[6 * 3] * (1 - ksi)*(eta)*(phi)+HexPnt[7 * 3] * (ksi)*(eta)*(phi);
		gc[1] = HexPnt[0 * 3 + 1] * (1 - ksi)*(1 - eta)*(1 - phi) + HexPnt[1 * 3 + 1] * (ksi)*(1 - eta)*(1 - phi) +
			HexPnt[2 * 3 + 1] * (1 - ksi)*(eta)*(1 - phi) + HexPnt[3 * 3 + 1] * (ksi)*(eta)*(1 - phi) +
			HexPnt[4 * 3 + 1] * (1 - ksi)*(1 - eta)*(phi)+HexPnt[5 * 3 + 1] * (ksi)*(1 - eta)*(phi)+
			HexPnt[6 * 3 + 1] * (1 - ksi)*(eta)*(phi)+HexPnt[7 * 3 + 1] * (ksi)*(eta)*(phi);
		gc[2] = HexPnt[0 * 3 + 2] * (1 - ksi)*(1 - eta)*(1 - phi) + HexPnt[1 * 3 + 2] * (ksi)*(1 - eta)*(1 - phi) +
			HexPnt[2 * 3 + 2] * (1 - ksi)*(eta)*(1 - phi) + HexPnt[3 * 3 + 2] * (ksi)*(eta)*(1 - phi) +
			HexPnt[4 * 3 + 2] * (1 - ksi)*(1 - eta)*(phi)+HexPnt[5 * 3 + 2] * (ksi)*(1 - eta)*(phi)+
			HexPnt[6 * 3 + 2] * (1 - ksi)*(eta)*(phi)+HexPnt[7 * 3 + 2] * (ksi)*(eta)*(phi);

	}

	bool calc_local_coordinates_in_hexahedron(double *R, double *HexPnt, double *lc)
	{
		int i;
		double EpsForFindLocalCoord = 1e-2;
		double CoeffForDiv = 1.2;
		int MaxDeep = 10;

		int p, t, m, deep, crd[3][2] = { { 1, 2 }, { 0, 2 }, { 0, 1 } };
		double LocalCoord[2][3], CentGlob[3], CentLoc[3], DopLoc[3];
		double dist, disc, h[3];
		const double ods = 0.16666666666666666;
		const double hds = 0.33333333333333333;

		double bb[6];
		bb[0] = HexPnt[0];
		bb[1] = HexPnt[0];
		bb[2] = HexPnt[1];
		bb[3] = HexPnt[1];
		bb[4] = HexPnt[2];
		bb[5] = HexPnt[2];

		for (i = 1; i < 6; i++)
		{
			if (bb[0] > HexPnt[i * 3])	bb[0] = HexPnt[i * 3];
			if (bb[1] < HexPnt[i * 3])	bb[1] = HexPnt[i * 3];
			if (bb[2] > HexPnt[i * 3 + 1])	bb[2] = HexPnt[i * 3 + 1];
			if (bb[3] < HexPnt[i * 3 + 1])	bb[3] = HexPnt[i * 3 + 1];
			if (bb[4] > HexPnt[i * 3 + 2])	bb[4] = HexPnt[i * 3 + 2];
			if (bb[5] < HexPnt[i * 3 + 2])	bb[5] = HexPnt[i * 3 + 2];
		}

		if (bb[0] > R[0])	return false;
		if (bb[1] < R[0])	return false;
		if (bb[2] > R[1])	return false;
		if (bb[3] < R[1])	return false;
		if (bb[4] > R[2])	return false;
		if (bb[5] < R[2])	return false;


		deep = 0;
		LocalCoord[0][0] = 0.0; LocalCoord[0][1] = 0.0; LocalCoord[0][2] = 0.0;
		LocalCoord[1][0] = 1.0; LocalCoord[1][1] = 1.0; LocalCoord[1][2] = 1.0;
		CentLoc[0] = 0.5; CentLoc[1] = 0.5; CentLoc[2] = 0.5;
		GetGlobalCoordinates(HexPnt, CentLoc, CentGlob);

		disc = sqrt(
			(R[0] - CentGlob[0])*(R[0] - CentGlob[0]) +
			(R[1] - CentGlob[1])*(R[1] - CentGlob[1]) +
			(R[2] - CentGlob[2])*(R[2] - CentGlob[2])
			);

		h[0] = LocalCoord[1][0] - LocalCoord[0][0];
		h[1] = LocalCoord[1][1] - LocalCoord[0][1];
		h[2] = LocalCoord[1][2] - LocalCoord[0][2];

		do
		{
			if (disc < EpsForFindLocalCoord)	break;

			for (m = 0; m < 3; m++)
			{

				CentLoc[m] = LocalCoord[0][m] + ods*h[m];

				GetGlobalCoordinates(HexPnt, CentLoc, CentGlob);
				dist = sqrt(
					(R[0] - CentGlob[0])*(R[0] - CentGlob[0]) +
					(R[1] - CentGlob[1])*(R[1] - CentGlob[1]) +
					(R[2] - CentGlob[2])*(R[2] - CentGlob[2])
					);

				if (disc < dist)
				{
					CentLoc[m] = LocalCoord[1][m] - ods*h[m];

					GetGlobalCoordinates(HexPnt, CentLoc, CentGlob);
					dist = sqrt(
						(R[0] - CentGlob[0])*(R[0] - CentGlob[0]) +
						(R[1] - CentGlob[1])*(R[1] - CentGlob[1]) +
						(R[2] - CentGlob[2])*(R[2] - CentGlob[2])
						);

					if (dist < disc)
					{
						disc = dist;
						LocalCoord[0][m] = LocalCoord[1][m] - hds*h[m];
					}
					else
					{
						LocalCoord[0][m] = LocalCoord[0][m] + hds*h[m];
						LocalCoord[1][m] = LocalCoord[1][m] - hds*h[m];
					}
				}
				else
				{
					disc = dist;
					LocalCoord[1][m] = LocalCoord[0][m] + hds*h[m];
				}
				h[m] = LocalCoord[1][m] - LocalCoord[0][m];
				CentLoc[m] = 0.5*(LocalCoord[0][m] + LocalCoord[1][m]);
			}

			deep++;
		}

		while (deep < MaxDeep);

		// Допоиск
		if (deep == MaxDeep)
		{
			DopLoc[0] = CentLoc[0];
			DopLoc[1] = CentLoc[1];
			DopLoc[2] = CentLoc[2];
			do
			{
				t = 0;
				for (m = 0; m < 3; m++)
				{
					do
					{
						p = 0;
						DopLoc[m] = CentLoc[m] - h[m];
						if (DopLoc[m] > 0){
							GetGlobalCoordinates(HexPnt, DopLoc, CentGlob);
							dist = sqrt(
								(R[0] - CentGlob[0])*(R[0] - CentGlob[0]) +
								(R[1] - CentGlob[1])*(R[1] - CentGlob[1]) +
								(R[2] - CentGlob[2])*(R[2] - CentGlob[2])
								);
							if (dist < disc){
								disc = dist;
								CentLoc[m] = DopLoc[m];
								t = p = 1;
							}
						}
					} while (p);
					do
					{
						p = 0;
						DopLoc[m] = CentLoc[m] + h[m];
						if (DopLoc[m] < 1){
							GetGlobalCoordinates(HexPnt, DopLoc, CentGlob);
							dist = sqrt(
								(R[0] - CentGlob[0])*(R[0] - CentGlob[0]) +
								(R[1] - CentGlob[1])*(R[1] - CentGlob[1]) +
								(R[2] - CentGlob[2])*(R[2] - CentGlob[2])
								);
							if (dist < disc){
								disc = dist;
								CentLoc[m] = DopLoc[m];
								t = p = 1;
							}
						}
					} while (p);
				}
			} while (t);
		}

		GetGlobalCoordinates(HexPnt, CentLoc, CentGlob);

		if (disc > EpsForFindLocalCoord)
			return false;

		if (CentLoc[0] < 0.0) return false;
		if (CentLoc[1] < 0.0) return false;
		if (CentLoc[2] < 0.0) return false;
		if (CentLoc[0] > 1.0) return false;
		if (CentLoc[1] > 1.0) return false;
		if (CentLoc[2] > 1.0) return false;

		lc[0] = CentLoc[0];
		lc[1] = CentLoc[1];
		lc[2] = CentLoc[2];

		return true;
	}

	int open_log(char *file_name)
	{
		if (open_file_w(file_name, &log_file) != 0)
		{
			printf("Error : could not open file '%s'\n", file_name);
			return 1;
		}

		return 0;
	}

	void write_to_log(char *str)
	{
		printf("%s", str);
		fprintf(log_file, "%s", str);
		fflush(log_file);
	}
	void write_to_log(const char *str)
	{
		printf("%s", str);
		fprintf(log_file, "%s", str);
		fflush(log_file);
	}
