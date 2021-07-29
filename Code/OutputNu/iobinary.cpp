#include "stdafx.h"
#include "iobinary.h"

ifstream& operator>(ifstream& file, short &id)
{
	file.read((char*)&id, sizeof(short));
	return file;
}

ofstream& operator<(ofstream& file, const short &id)
{
	file.write((char*)&id, sizeof(short));
	return file;
}

ifstream& operator>(ifstream& file, int &id)
{
	file.read((char*)&id, sizeof(int));
	return file;
}

ofstream& operator<(ofstream& file, const int &id)
{
	file.write((char*)&id, sizeof(int));
	return file;
}

ifstream& operator>(ifstream& file, long &id)
{
	file.read((char*)&id, sizeof(long));
	return file;
}

ofstream& operator<(ofstream& file, const long &id)
{
	file.write((char*)&id, sizeof(long));
	return file;
}

ifstream& operator>(ifstream& file, double &id)
{
	file.read((char*)&id, sizeof(double));
	return file;
}

ofstream& operator<(ofstream& file, const double &id)
{
	file.write((char*)&id, sizeof(double));
	return file;
}

ifstream& operator>(ifstream& file, float &id)
{
	file.read((char*)&id, sizeof(float));
	return file;
}

ofstream& operator<(ofstream& file, const float &id)
{
	file.write((char*)&id, sizeof(float));
	return file;
}
