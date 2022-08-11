// BasicOperations.cpp header file
//
// Phd. Kiselev D.S.
// Novosibirsk State Technical University
// 20, Karl Marx Avenue, Novosibirsk, Russia 

#ifndef __BasicOperations__H__
#define __BasicOperations__H__
#include <set>
#include <map>

extern FILE *log_file;
using namespace std;
#define _PI_ 3.14159265358979

	int open_file_w(char *file_name, FILE **file_stream);

	int open_log(char *file_name);

	void write_to_log(char *str);
	void write_to_log(const char *str);
#endif