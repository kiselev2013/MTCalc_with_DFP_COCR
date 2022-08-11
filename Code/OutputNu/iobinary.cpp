/**
 * GENERAL REMARKS
 *
 *  This code is freely available under the following conditions:
 *
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer.
 *
 *              MTCalc_with_DFP_COCR
 *  Functions for binary input-output
 *
 *  Written by Prof. Marina G. Persova
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  mpersova@mail.ru
 *  Version 1.5 December 17, 2020
*/
#include "stdafx.h"
#include "iobinary.h"
//-----------------------------------------------------------
// Binary reading for element of type short int
//-----------------------------------------------------------
ifstream& operator>(ifstream& file, short &id)
{
	file.read((char*)&id, sizeof(short));
	return file;
}
//-----------------------------------------------------------
// Binary writing for element of type short int
//-----------------------------------------------------------
ofstream& operator<(ofstream& file, const short &id)
{
	file.write((char*)&id, sizeof(short));
	return file;
}
//-----------------------------------------------------------
// Binary reading for element of type int
//-----------------------------------------------------------
ifstream& operator>(ifstream& file, int &id)
{
	file.read((char*)&id, sizeof(int));
	return file;
}
//-----------------------------------------------------------
// Binary writing for element of type int
//-----------------------------------------------------------
ofstream& operator<(ofstream& file, const int &id)
{
	file.write((char*)&id, sizeof(int));
	return file;
}
//-----------------------------------------------------------
// Binary reading for element of type long
//-----------------------------------------------------------
ifstream& operator>(ifstream& file, long &id)
{
	file.read((char*)&id, sizeof(long));
	return file;
}
//-----------------------------------------------------------
// Binary writing for element of type long
//-----------------------------------------------------------
ofstream& operator<(ofstream& file, const long &id)
{
	file.write((char*)&id, sizeof(long));
	return file;
}
//-----------------------------------------------------------
// Binary reading for element of type double
//-----------------------------------------------------------
ifstream& operator>(ifstream& file, double &id)
{
	file.read((char*)&id, sizeof(double));
	return file;
}
//-----------------------------------------------------------
// Binary writing for element of type double
//-----------------------------------------------------------
ofstream& operator<(ofstream& file, const double &id)
{
	file.write((char*)&id, sizeof(double));
	return file;
}
//-----------------------------------------------------------
// Binary reading for element of type float
//-----------------------------------------------------------
ifstream& operator>(ifstream& file, float &id)
{
	file.read((char*)&id, sizeof(float));
	return file;
}
//-----------------------------------------------------------
// Binary writing for element of type float
//-----------------------------------------------------------
ofstream& operator<(ofstream& file, const float &id)
{
	file.write((char*)&id, sizeof(float));
	return file;
}
