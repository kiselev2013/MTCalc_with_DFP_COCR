#include "stdafx.h"
#include "divgrad_1d.h"

ofstream logfile;

int main()
{
	logfile.open("LogMTZ1D");

	Divgrad_1d divgrad1d;

	divgrad1d.Solve_1d_Problem_for_3d_task();

	logfile.close();
	logfile.clear();

	return 0;
}