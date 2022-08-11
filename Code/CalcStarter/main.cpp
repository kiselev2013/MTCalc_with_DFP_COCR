/**
 *  This code is freely available under the following conditions:
 *
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer.
 *
 *
 * The programm runs main calculation procedures sequense
 * Contains functions for preparing calculation directories and 
 * starting calculation programs
 *
 *  Written by Ph.D. Dmitry S. Kiselev
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  harlequin_00@mail.ru
 *  Version 1.1 July, 2021
*/

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <Windows.h>
#include <direct.h>

using namespace std;

char currentBigLetter[] = {'Y', 'X'};
char currentSmallLetter[] = { 'y', 'x' };

// Calculation settings class
class Settings
{
public:
	double steps[3]; // Mesh initial steps
	double sparse[3]; // Sparce coefficient
	double farBound[3]; // Distance from receivers bounding box+gap to the boundary of the calculation domain
	double gapFromReceivers[3]; // Gap between receivers bounding box and sparced mesh area
	double eps = 0.1; // Gate for determining if coordinate lines are equal
	int threadsCount = 1; // Threads count for parallel calculation
	double SLAESolutionEps = 1e-4; // relative residual for determining if SLAE solving is finished
	int maxIterationsCount = 10000; // Maximum iterations when solving SLAE

	Settings() { }
};

FILE *log_file = NULL;

// Write message to log
void write_to_log(char *str)
{
	printf("%s", str);
	if (log_file != NULL)
	{
		fprintf(log_file, "%s", str);
		fflush(log_file);
	}
}

// Write message to log
void write_to_log(const char *str)
{
	printf("%s", str);
	if (log_file != NULL)
	{
		fprintf(log_file, "%s", str);
		fflush(log_file);
	}
}

// Open file for reading
int open_file_r(char *file_name, FILE **file_stream)
{
	char buf[2048];

	if (!((*file_stream) = fopen(file_name, "r")))
	{
		sprintf(buf, "Error : Could not read file '%s'\n", file_name);
		write_to_log(buf);
		return 1;
	}

	return 0;
}

// Open file for writing 
int open_file_w(char *file_name, FILE **file_stream)
{
	char buf[2048];

	if (!((*file_stream) = fopen(file_name, "w")))
	{
		sprintf(buf, "Error : Could not write file '%s'\n", file_name);
		write_to_log(buf);
		return 1;
	}

	return 0;
}

// Open log file
int open_log(char *file_name)
{
	if (open_file_w(file_name, &log_file) != 0)
	{
		printf("Error : could not open file '%s'\n", file_name);
		return 1;
	}

	return 0;
}

// Run Windows executable file
int ExecuteExe(const char *cmdline, char *workdir)
{
	int retCode;
	STARTUPINFOA si;
	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);

	PROCESS_INFORMATION pi;
	ZeroMemory(&pi, sizeof(pi));

	CreateProcessA(NULL, (LPSTR)cmdline, NULL, NULL, FALSE, CREATE_NO_WINDOW, NULL, workdir, &si, &pi);
	WaitForSingleObject(pi.hProcess, INFINITE);
	GetExitCodeProcess(pi.hProcess, (LPDWORD)&retCode);

	CloseHandle(pi.hProcess);
	CloseHandle(pi.hThread);
	return retCode;
}

// Read application settings file
int ReadSettings(char *file_name, Settings &settings)
{
	FILE *file_in = NULL;
	if (open_file_r(file_name, &file_in) != 0) return 1;
	char c;

	fscanf(file_in, "%lf%lf", settings.steps, settings.steps + 2); settings.steps[1] = settings.steps[0];
	fscanf(file_in, "%lf%lf", settings.sparse, settings.sparse + 2); settings.sparse[1] = settings.sparse[0];
	fscanf(file_in, "%lf%lf", settings.gapFromReceivers, settings.gapFromReceivers + 2); settings.gapFromReceivers[1] = settings.gapFromReceivers[0];
	fscanf(file_in, "%lf%lf", settings.farBound, settings.farBound + 2); settings.farBound[1] = settings.farBound[0];
	fscanf(file_in, "%lf", &settings.SLAESolutionEps);
	fscanf(file_in, "%d", &settings.maxIterationsCount);

	fclose(file_in);
	return 0;
}

// Read frequencies
int ReadFrequencies(char *file_name, vector<double> &frequencies)
{
	FILE *file_in = NULL;
	if (open_file_r(file_name, &file_in) != 0) return 1;

	int frequenciesCount;
	if (fscanf(file_in, "%d", &frequenciesCount) != 1)
	{
		fclose(file_in);
		return 1;
	}

	frequencies.resize(frequenciesCount);
	for (int frequencyIndex = 0; frequencyIndex < frequenciesCount; frequencyIndex++)
	{
		if (fscanf(file_in, "%lf", &frequencies[frequencyIndex]) != 1)
		{
			fclose(file_in);
			return 1;
		}
	}

	fclose(file_in);
	return 0;
}

// Write settings file
int WriteSettings(char *file_name, Settings &settings, int currentDirection, double frequency)
{
	FILE *file_out = NULL;
	if (open_file_w(file_name, &file_out) != 0) return 1;
	char c;

	fprintf(file_out, "%.13e\t%.13e\n", settings.steps[0], settings.steps[2]);
	fprintf(file_out, "%.13e\t%.13e\n", settings.sparse[0], settings.sparse[2]);
	fprintf(file_out, "%.13e\t%.13e\n", settings.gapFromReceivers[0], settings.gapFromReceivers[2]);
	fprintf(file_out, "%.13e\t%.13e\n", settings.farBound[0], settings.farBound[2]);
	fprintf(file_out, "%c\n", currentBigLetter[currentDirection]);
	fprintf(file_out, "%.13e\n", frequency);
	fprintf(file_out, "%d\n", 1);
	fprintf(file_out, "%.13e\n", settings.SLAESolutionEps);
	fprintf(file_out, "%d\n", settings.maxIterationsCount);

	fclose(file_out);
	return 0;

}

// Write frequencies file
int WriteFreq(char *file_name, vector<double> &frequencies)
{
	FILE *file_out = NULL;
	if (open_file_w(file_name, &file_out) != 0) return 1;

	for (int frequencyIndex = 0; frequencyIndex < frequencies.size(); frequencyIndex++)
		fprintf(file_out, "%.13e\n", frequencies[frequencyIndex]);

	fclose(file_out);
	return 0;
}

// Write frequencies file (another format)
int WriteFrec(char *file_name, vector<double> &frequencies)
{
	FILE *file_out = NULL;
	if (open_file_w(file_name, &file_out) != 0) return 1;

	fprintf(file_out, "%d\n", frequencies.size());
	for (int frequencyIndex = frequencies.size()-1; frequencyIndex >= 0; frequencyIndex--)
		fprintf(file_out, "%.13e\n", frequencies[frequencyIndex]);

	fclose(file_out);
	return 0;
}

// Prepare directory for frequency, current direction
int PrepareDirectoryNu(Settings &settings, int frequencyIndex, double frequency, int currentDirection)
{
	char path[2048], fileName[2048];
	sprintf(path, "Calculation/Nu%d_%c", frequencyIndex + 1, currentBigLetter[currentDirection]);
	_mkdir(path);

	sprintf(fileName, "%s/objects", path); CopyFileA("objects", fileName, false);
	sprintf(fileName, "%s/z_sig_2d", path); CopyFileA("z_sig_2d", fileName, false);
	sprintf(fileName, "%s/xyzVectorB", path); CopyFileA("xyzVectorB", fileName, false);

	sprintf(fileName, "%s/RegularMeshBuilderSettings.cfg", path);
	if (WriteSettings(fileName, settings, currentDirection, frequency) != 0)
		return 1;

	return 0;
}

// Prepare directory for frequency calculation
int PrepareDirectory(Settings &settings, vector<double> &frequencies)
{
	_mkdir("Calculation");
	_mkdir("Results");

	for (int frequencyIndex = 0; frequencyIndex < frequencies.size(); frequencyIndex++)
	{
		if (PrepareDirectoryNu(settings, frequencyIndex, frequencies[frequencyIndex], 0) != 0) return 1;
		if (PrepareDirectoryNu(settings, frequencyIndex, frequencies[frequencyIndex], 1) != 0) return 1;
	}

	return 0;
}

// Copying result files to common 'Results' folder
int CopyResultFiles(int frequencyIndex, int currentDirection)
{
	char path[2048], fileName[2048], sourceFileName[2048], destFileName[2048], buf[2048];

	sprintf(path, "Calculation/Nu%d_%c", frequencyIndex + 1, currentBigLetter[currentDirection]);

	sprintf(fileName, "ex_s"); sprintf(sourceFileName, "%s/%s", path, fileName); sprintf(destFileName, "Results/%s%c_s0f%d", fileName, currentBigLetter[currentDirection], frequencyIndex); CopyFileA(sourceFileName, destFileName, false);
	sprintf(fileName, "ex_c"); sprintf(sourceFileName, "%s/%s", path, fileName); sprintf(destFileName, "Results/%s%c_s0f%d", fileName, currentBigLetter[currentDirection], frequencyIndex); CopyFileA(sourceFileName, destFileName, false);
	sprintf(fileName, "ey_s"); sprintf(sourceFileName, "%s/%s", path, fileName); sprintf(destFileName, "Results/%s%c_s0f%d", fileName, currentBigLetter[currentDirection], frequencyIndex); CopyFileA(sourceFileName, destFileName, false);
	sprintf(fileName, "ey_c"); sprintf(sourceFileName, "%s/%s", path, fileName); sprintf(destFileName, "Results/%s%c_s0f%d", fileName, currentBigLetter[currentDirection], frequencyIndex); CopyFileA(sourceFileName, destFileName, false);
	sprintf(fileName, "ez_s"); sprintf(sourceFileName, "%s/%s", path, fileName); sprintf(destFileName, "Results/%s%c_s0f%d", fileName, currentBigLetter[currentDirection], frequencyIndex); CopyFileA(sourceFileName, destFileName, false);
	sprintf(fileName, "ez_c"); sprintf(sourceFileName, "%s/%s", path, fileName); sprintf(destFileName, "Results/%s%c_s0f%d", fileName, currentBigLetter[currentDirection], frequencyIndex); CopyFileA(sourceFileName, destFileName, false);
	sprintf(fileName, "hx_s"); sprintf(sourceFileName, "%s/%s", path, fileName); sprintf(destFileName, "Results/%s%c_s0f%d", fileName, currentBigLetter[currentDirection], frequencyIndex); CopyFileA(sourceFileName, destFileName, false);
	sprintf(fileName, "hx_c"); sprintf(sourceFileName, "%s/%s", path, fileName); sprintf(destFileName, "Results/%s%c_s0f%d", fileName, currentBigLetter[currentDirection], frequencyIndex); CopyFileA(sourceFileName, destFileName, false);
	sprintf(fileName, "hy_s"); sprintf(sourceFileName, "%s/%s", path, fileName); sprintf(destFileName, "Results/%s%c_s0f%d", fileName, currentBigLetter[currentDirection], frequencyIndex); CopyFileA(sourceFileName, destFileName, false);
	sprintf(fileName, "hy_c"); sprintf(sourceFileName, "%s/%s", path, fileName); sprintf(destFileName, "Results/%s%c_s0f%d", fileName, currentBigLetter[currentDirection], frequencyIndex); CopyFileA(sourceFileName, destFileName, false);
	sprintf(fileName, "hz_s"); sprintf(sourceFileName, "%s/%s", path, fileName); sprintf(destFileName, "Results/%s%c_s0f%d", fileName, currentBigLetter[currentDirection], frequencyIndex); CopyFileA(sourceFileName, destFileName, false);
	sprintf(fileName, "hz_c"); sprintf(sourceFileName, "%s/%s", path, fileName); sprintf(destFileName, "Results/%s%c_s0f%d", fileName, currentBigLetter[currentDirection], frequencyIndex); CopyFileA(sourceFileName, destFileName, false);
	sprintf(fileName,  "e2d"); sprintf(sourceFileName, "%s/%s", path, fileName); sprintf(destFileName, "Results/%s%c_s0f%d", fileName, currentBigLetter[currentDirection], frequencyIndex); CopyFileA(sourceFileName, destFileName, false);
	sprintf(fileName,  "b2d"); sprintf(sourceFileName, "%s/%s", path, fileName); sprintf(destFileName, "Results/%s%c_s0f%d", fileName, currentBigLetter[currentDirection], frequencyIndex); CopyFileA(sourceFileName, destFileName, false);

	if (frequencyIndex == 0 && currentDirection == 1)
	{
		sprintf(fileName, "pointres"); sprintf(sourceFileName, "%s/%s", path, fileName); sprintf(destFileName, "Results/%s", fileName); CopyFileA(sourceFileName, destFileName, false);
	}

	return 0;
}

// Run calculation for frequency, current direction
int RunCalculation(int frequencyIndex, int currentDirection)
{
	char path[2048], buf[2048];
	
	sprintf(path, "Calculation/Nu%d_%c", frequencyIndex + 1, currentBigLetter[currentDirection]);
	sprintf(buf, "Starting calculation in %s\n", path);
	write_to_log(buf);

	if (ExecuteExe("Programs/RegularMeshBuilder.exe", path) != 0) return 1;
	if (ExecuteExe("Programs/Harm1D.exe", path) != 0) return 1;
	if (ExecuteExe("Programs/BuildMatrix.exe", path) != 0) return 1;
	if (ExecuteExe("Programs/COCR_FP.exe", path) != 0) return 1;
	if (ExecuteExe("Programs/OutputNu.exe", path) != 0) return 1;

	CopyResultFiles(frequencyIndex, currentDirection);

	return 0;
}

// Run calculations for all frequencies
int RunCalculations(vector<double> &frequencies)
{
	for (int frequencyIndex = 0; frequencyIndex < frequencies.size(); frequencyIndex++)
	{
		if (RunCalculation(frequencyIndex, 1) != 0) return 1;
		if (RunCalculation(frequencyIndex, 0) != 0) return 1;
	}

	if (WriteFreq("Results/freq", frequencies) != 0) return 1;
	if (WriteFrec("Results/frec", frequencies) != 0) return 1;

	if (ExecuteExe("Programs/ImpedanceZToEdsAll.exe", "Results") != 0) return 1;

	return 0;
}

// Calling common procedures
int MainProcedure()
{
	Settings settings;
	vector<double> frequencies;

	if (ReadSettings("settings.cfg", settings) != 0)
		return 1;

	if (ReadFrequencies("frequencies", frequencies) != 0)
		return 1;

	if (PrepareDirectory(settings, frequencies) != 0)
		return 1;

	if (RunCalculations(frequencies) != 0)
		return 1;

	return 0;
}

int main()
{
	char buf[2048];

	if (open_log("CalcStarter.log") != 0)
	{
		printf("Cound not open CalcStarter.log\n");
		return 1;
	}

	int status = MainProcedure();
	sprintf(buf, "Status = %d\n", status);
	write_to_log(buf);

	fprintf(log_file, "All done\n");
	fclose(log_file);
	return status;
}
