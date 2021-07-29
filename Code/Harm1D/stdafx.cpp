#include "stdafx.h"

extern ofstream logfile;

void Memory_allocation_error(const char *var, const char *func)
{
	string str;
	str = "MEMORY ALLOCATION ERROR for variable ";
	str = str + '\"' + var + '\"' + " in function " + '\"' + func + "\"\n";
	cerr    << str << flush;
	cout << str << flush;
	throw logic_error(str);
}
//-----------------------------------------------------------
void Cannot_open_file(const char *fname, const char *func)
{
	string str;
	str = "CANNOT OPEN FILE ";
	str = str + '\"' + fname + '\"' + " in function " + '\"' + func + "\"\n";
	cerr    << str << flush;
	cout << str << flush;
	throw logic_error(str);
}
//-----------------------------------------------------------
void Cannot_open_file_but_continue(const char *fname, const char *func)
{
	string str;
	str = "Cannot open file ";
	str = str + '\"' + fname + '\"' + " in function " + '\"' + func + "\"\n";
	cerr    << str << flush;
}