// To compile, use Microsoft Visual Studio C++ compiler
//
// This module contains some utils and logging function
//
// Phd. Persova M.G.
// Novosibirsk State Technical University
// 20, Karl Marx Avenue, Novosibirsk, Russia 

#include "stdafx.h"
#include "utils.h"

// Writing error message
void WriteMessage(char *str)
{
	cout<<str;
	logfile<<str;
}

// Writing error message
void WriteMessage(stringstream &str)
{
	cout<<str.str();
	logfile<<str.str();
}

// Checking if a file exists
bool IsFileExist(char *fname)
{
	bool flag;
	ifstream inf;

	flag=false;
	inf.open(fname);
	if(inf)
	{
		flag=true;
		inf.close();
	}
	inf.clear();

	return flag;
}

// Closing programm with error
void CloseProgramm(int err_code)
{
	stringbuffer<<"Closing programm with err_code "<<err_code<<endl;
	exit(err_code);
}

// Checking if a error occured while reading file
void StopIfErrorReturn(int err_code,char *FuncName)
{
	if(err_code)
	{
		stringbuffer<<"Function "<<FuncName<<" returned "<<err_code<<endl;
		WriteMessage(stringbuffer);
		CloseProgramm(err_code);
	}
}
