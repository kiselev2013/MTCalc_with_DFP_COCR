// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include <stdio.h>
#include <tchar.h>


// TODO: reference additional headers your program requires here
#include <math.h>
#include <omp.h>
#include <time.h>


#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stack>
#include <complex>
#include <iomanip>
using namespace std;

#include "in_out.h"
#include "For_Solvers.h"


#define MU_0  1.25663706143591729e-6
#define PI    3.1415926535897932384626433832795
#define DPR_0 8.84194128288307421e-12

#define _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_DEPRECATE


__forceinline ostream& operator < (ostream& file,const double& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,double&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const float& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,float&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const int& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,int&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const bool& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,bool&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const char& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,char&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const unsigned char& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,unsigned char&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const short& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,short&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const unsigned int & data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,unsigned int &  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}


