#include "stdafx.h"
#include "open_close.h"

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
