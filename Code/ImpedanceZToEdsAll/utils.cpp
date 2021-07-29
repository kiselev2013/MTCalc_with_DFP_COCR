#include "stdafx.h"
#include "utils.h"

void WriteMessage(char *str)
{
	cout<<str;
	logfile<<str;
}

void WriteMessage(stringstream &str)
{
	cout<<str.str();
	logfile<<str.str();
}

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

void CloseProgramm(int err_code)
{
	stringbuffer<<"Closing programm with err_code "<<err_code<<endl;
	exit(err_code);
}

void StopIfErrorReturn(int err_code,char *FuncName)
{
	if(err_code)
	{
		stringbuffer<<"Function "<<FuncName<<" returned "<<err_code<<endl;
		WriteMessage(stringbuffer);
		CloseProgramm(err_code);
	}
}

int CreateProcessForEXENoWait(char *cmdline, char *workdir)
{
	int retp;
	STARTUPINFOA si;
	PROCESS_INFORMATION pi;
	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);
	ZeroMemory(&pi, sizeof(pi));
	cout<<"Start "<<cmdline;
	if(workdir){cout<<" in "<<workdir;}
	cout<<endl;
	// Start the child process. 
	retp=CreateProcessA(NULL,  // No module name (use command line). 
		(LPSTR)(const char*)cmdline,// Command line. 
		NULL,				// Process handle not inheritable. 
		NULL,				// Thread handle not inheritable. 
		FALSE,				// Set handle inheritance to FALSE. 
		CREATE_NO_WINDOW,	// No creation flags. 
		NULL,				// Use parent's environment block. 
		workdir,				// Use parent's starting directory. 
		&si,				// Pointer to STARTUPINFO structure.
		&pi);				// Pointer to PROCESS_INFORMATION structure.
	return 0;
}
