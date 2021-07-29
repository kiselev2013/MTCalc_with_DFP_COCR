#include "stdafx.h"
#include "linalg.h"
#include "retcode.h"
//------------------------------------------------------------------------
// Скалярное произведение векторов
//------------------------------------------------------------------------
complex<double> dot(complex<double> *a, complex<double> *b, int n)
{
	complex<double> s = 0;
	int i;
	int m = n%4;

	for (i=0; i<m; i++)
	{
		s += a[i]*b[i];
	}

	for (i=m; i<n; i+=4)
	{
		s += a[i]*b[i] + a[i+1]*b[i+1] + a[i+2]*b[i+2] + a[i+3]*b[i+3];
	}

	return s;
}
void ger(complex<double> *a, complex<double> *x, complex<double> *y, int n, int m, int lda)
{
	

	int i, j;
	int adr;
	complex<double> temp;
	int k = m%4;

	#pragma omp parallel for private(adr, j, temp)
	for (i=0; i<n; i++)
	{
		adr = (lda-n+i)*lda + lda-n;
		temp = x[i*lda];

		for (j=0; j<k; j++)
		{
			a[adr + j] -= y[j]*temp;
		}

		for (j=k; j<m; j+=4)
		{
			a[adr + j]     -= y[j]*temp;
			a[adr + j + 1] -= y[j+1]*temp;
			a[adr + j + 2] -= y[j+2]*temp;
			a[adr + j + 3] -= y[j+3]*temp;
		}
	}
}
//------------------------------------------------------------------------
// LU-разложение
//------------------------------------------------------------------------
int getf2(complex<double> *a, int m, int n, int lda)
{
	int i, j;
	int dim;
	int nstr;
	int adr;
	complex<double> temp;

	dim = m - 1;
	if (n < dim)
		dim = n;

	for (i=0; i<dim; i++)
	{
		nstr = lda - m + i;
		adr = nstr*lda + nstr;
		temp = a[adr];

		if (abs(temp) < 1e-30)
			return -1;

		temp = 1.0/temp;

		for (j=1; j<m-i; j++)
			a[adr + j*lda] *= temp;

		if (i < n-1)
			ger(a, &a[(nstr + 1)*lda + nstr], &a[nstr*lda + nstr + 1], m - i - 1, n - i - 1, lda);
	}	

	return 0;
}
//------------------------------------------------------------------------
// Решение СЛАУ с нижнетреугольной матрицей, на диагонали 1
//------------------------------------------------------------------------
void SolveL1(complex<double> *a, complex<double> *b, int n)
{
	int i;
	complex<double> s;

	for (i=0; i<n; i++)
	{
		s = dot(&a[i*n], b, i);
		b[i] -= s;
	}
}
//------------------------------------------------------------------------
// Решение СЛАУ с верхнетреугольной матрицей
//------------------------------------------------------------------------
int SolveU(complex<double> *a, complex<double> *b, int n)
{
	int i,ret;
	complex<double> s;
	ret=RETCODE_OK;
	for (i=n-1; i>=0; i--)
	{
		s = dot(&a[i*n+i+1], &b[i+1], n-1-i);
		if(abs(a[i*n + i])<1e-20)
		{
			ret=RETCODE_DEVBYZERO;
			break;
		}
		b[i] = (b[i] - s)/a[i*n + i];
	}
	return ret;
}
//------------------------------------------------------------------------
// LU-разложение
//------------------------------------------------------------------------
int getrf(complex<double> *a, int n,const int b,int blocksize,complex<double> *L,complex<double> *f,complex<double> *aa)
{
	int flag1 = 1;
	int flag2 = 1;
	int i;

	if (b <= 1 || b >= n)
	{
		return getf2(a, n, n, n);
	} 
	else
	{
		for (i=0; i<n; i+=b)
		{
			blocksize = n-i;
			if (b < blocksize)
				blocksize = b;

			if (getf2(a,  n-i, blocksize, n)!=0)
				break;



			pre_trsm(a, L, f, blocksize, n,  i);
			trsm(L, f, blocksize, n - i - blocksize);
			post_trsm(a, f, blocksize, n, i);



			pre_MM(a, aa, blocksize, i, n);

	 		if (blocksize == b)
	 		{
				MultMVvectorizeBlock(aa, f, a, n-i-blocksize, blocksize,  n);
	 		} 
	 		else
	 		{
				MultMVblock(aa, f, a, n-i-blocksize, blocksize,  n);
	 		}
		}
	}

	return 0;
}
//------------------------------------------------------------------------
void pre_MM(complex<double> *a, complex<double> *b, int n, int m, int lda)
{
	int i, j;
	int adr;
	int kol = lda - n - m;

	for (i=0; i<kol; i++)
	{
		adr = (m + n + i)*lda + m;

		for (j=0; j<n; j++)
		{
			b[i*n + j] = a[adr + j];
		}
	}
}
//------------------------------------------------------------------------
void pre_trsm(complex<double> *a, complex<double> *L, complex<double> *b, int n, int lda, int m)
{
/*
	Перезаписывает нижний треугольник подматрицы "L" матрицы "a" в отдельный
	массив. Матрица "a" размера [lda]x[lda]. Матрица "L" размера [n]x[n]
	начинается с элемента [m, m] матрицы "a".

	Перезаписывает подматрицу a[m+n:lda][m:m+n] в массив "b", где элементы
	хранятся в столбцовом порядке.
*/
	int i, j, k;
	int adr;
	int kol = lda - m - n;

	k = 0;
	for (i=0; i<n; i++)
	{
		adr = (m + i)*lda + m;
		for(j=0; j<i; j++)
		{
			L[k] = a[adr + j];
			k++;
		}
	}


	for (i=0; i<n; i++)
	{
		adr = (m + i)*lda + m + n;

		for (j=0; j<kol; j++)
			b[j*n+i] = a[adr + j];
	}
}
//------------------------------------------------------------------------
void post_trsm(complex<double> *a, complex<double>*b, int n, int lda, int m)
{
	int i, j;
	int adr;

	for (i=0; i<n; i++)
	{
		adr = (m + i)*lda + m + n;

		for (j=0; j<lda-m-n; j++)
			a[adr + j] = b[j*n+i];
	}
}
//------------------------------------------------------------------------
void trsm(complex<double> *L, complex<double> *b, int n, int m)
{
	int i, j, k;
	complex<double> s;

	#pragma omp parallel private(i, j, k, s)
	{
		k = 0;

		for (i=0; i<n; i++)
		{
			#pragma omp for schedule(static) nowait
			for (j=0; j<m; j++)
			{
				s = dot(&L[k], &b[j*n], i);
				b[j*n + i] -= s;
			}

			k += i;
		}
	}
}
//------------------------------------------------------------------------
void MultMVblock(complex<double> *a, complex<double> *b, complex<double> *c, int n, int m, int lda)
{
	int nb = n/m;
	int ib;

	int i, j, k;
	int adr;

	#pragma omp parallel private (i, j, k, adr, ib)
	{
		#pragma omp for
		for (ib=0; ib<nb; ib++)
		{
			for (i=ib*m; i<ib*m+m; i++)
			{
				for (j=0; j<n; j++)		
				{
					adr = (lda-n + i)*lda + (lda-n + j);

					for (k=0; k<m; k++)
						c[adr] -= a[i*m + k]*b[j*m + k];
				}
			}
		}

		#pragma omp for
		for (i=nb*m; i<n; i++)
		{
			for (j=0; j<n; j++)		
			{
				adr = (lda-n + i)*lda + (lda-n + j);

				for (k=0; k<m; k++)
					c[adr] -= a[i*m + k]*b[j*m + k];
			}
		}
	}
}
//------------------------------------------------------------------------
void MultMV(complex<double> *a, complex<double> *b, complex<double> *c, int n, int m, int lda)
{
	int i, j, k;
	int adr;

	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)		
		{
			adr = (lda-n + i)*lda + (lda-n + j);

			for (k=0; k<m; k++)
				c[adr] -= a[i*m + k]*b[j*m + k];
		}
	}
}
//------------------------------------------------------------------------
void MultMVvectorize(complex<double> *a, complex<double> *b, complex<double> *c, int n, int m, int lda)
{
	int i, j, k;
	int adr;

	#pragma omp parallel for private (j, k, adr)

	for (i=0; i<n; i++)
	{
		for (j=0; j<n; j++)		
		{
			adr = (lda-n + i)*lda + (lda-n + j);
			for (k=0; k<m; k+=4)
			{
				c[adr] -= a[i*m + k]*b[j*m + k];
				c[adr] -= a[i*m + k+1]*b[j*m + k+1];
				c[adr] -= a[i*m + k+2]*b[j*m + k+2];
				c[adr] -= a[i*m + k+3]*b[j*m + k+3];
			}
		}
	}
}
//------------------------------------------------------------------------
void MultMVvectorizeBlock(complex<double> *a, complex<double> *b, complex<double> *c, int n, int m, int lda)
{
	int i, j, k;
	int adr;
	int nb = n/m;
	int ib;

	#pragma omp parallel private (i, j, k, adr, ib)
	{
		#pragma omp for
		for (ib=0; ib<nb; ib++)
		{
			for (i=ib*m; i<ib*m+m; i++)
			{
				for (j=0; j<n; j++)		
				{
					adr = (lda-n + i)*lda + (lda-n + j);
					for (k=0; k<m; k+=4)
					{
						c[adr] -= a[i*m + k]*b[j*m + k];
						c[adr] -= a[i*m + k+1]*b[j*m + k+1];
						c[adr] -= a[i*m + k+2]*b[j*m + k+2];
						c[adr] -= a[i*m + k+3]*b[j*m + k+3];
					}
				}
			}
		}

		
		#pragma omp for
		for (i=nb*m; i<n; i++)
		{
			for (j=0; j<n; j++)		
			{
				adr = (lda-n + i)*lda + (lda-n + j);
				for (k=0; k<m; k+=4)
				{
					c[adr] -= a[i*m + k]*b[j*m + k];
					c[adr] -= a[i*m + k+1]*b[j*m + k+1];
					c[adr] -= a[i*m + k+2]*b[j*m + k+2];
					c[adr] -= a[i*m + k+3]*b[j*m + k+3];
				}
			}
		}
	}
}
//------------------------------------------------------------------------







