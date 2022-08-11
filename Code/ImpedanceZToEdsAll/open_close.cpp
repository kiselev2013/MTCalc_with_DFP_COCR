// To compile, use Microsoft Visual Studio C++ compiler
//
// This module contains functions-shells for standard file I/O functions.
//
// Phd. Persova M.G.
// Novosibirsk State Technical University
// 20, Karl Marx Avenue, Novosibirsk, Russia 

#include "stdafx.h"
#include "open_close.h"

// Opening file for reading
int OpenInputFile(ifstream &inf,char *fname,ios_base::open_mode mode)
{
	inf.open(fname,mode);
	if(inf)
	{
		return 0;
	}
	else
	{
		inf.clear();
		cout<<"Can't open file "<<fname<<endl;
		return 1;
	}
}

// Opening file for reading
int OpenOutputFile(ofstream &ofp,char *fname,ios_base::open_mode mode)
{
	ofp.open(fname,mode);
	if(ofp)
	{
		return 0;
	}
	else
	{
		ofp.clear();
		cout<<"Can't open file "<<fname<<endl;
		return 1;
	}
}

// Closing file after reading
int CloseInputFile(ifstream &inf,char *fname)
{
	if(inf.good())
	{
		inf.close();
		inf.clear();
		return 0;
	}
	else
	{
		inf.clear();
		cout<<"Error while reading file "<<fname<<endl;
		return 1;
	}
}

// Closing file after reading
int CloseOutputFile(ofstream &ofp,char *fname)
{
	if(ofp.good())
	{
		ofp.close();
		ofp.clear();
		return 0;
	}
	else
	{
		ofp.clear();
		cout<<"Error while writing file "<<fname<<endl;
		return 1;
	}
}
