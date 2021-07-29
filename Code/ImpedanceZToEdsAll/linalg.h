#pragma once
//------------------------------------------------------------------------
complex<double> dot(complex<double> *a, complex<double> *b, int n);
void ger(complex<double> *a, complex<double> *x, complex<double> *y, int n, int m, int lda);
int getf2(complex<double> *a, int m, int n, int lda);
void SolveL1(complex<double> *a, complex<double> *b, int n);
int SolveU(complex<double> *a, complex<double> *b, int n);
int getrf(complex<double> *a, int n,const int b,int blocksize,complex<double> *L,complex<double> *f,complex<double> *aa);
void pre_trsm(complex<double> *a, complex<double> *L, complex<double> *b, int n, int lda, int m);
void post_trsm(complex<double> *a, complex<double>*b, int n, int lda, int m);
void trsm(complex<double> *L, complex<double> *b, int n, int m);
void MultMV(complex<double> *a, complex<double> *b, complex<double> *c, int n, int m, int lda);
void MultMVvectorize(complex<double> *a, complex<double> *b, complex<double> *c, int n, int m, int lda);
void pre_MM(complex<double> *a, complex<double> *b, int n, int m, int lda);
void MultMVblock(complex<double> *a, complex<double> *b, complex<double> *c, int n, int m, int lda);
void MultMVvectorizeBlock(complex<double> *a, complex<double> *b, complex<double> *c, int n, int m, int lda);