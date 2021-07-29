#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>

#define PI 3.1415926535897932

using namespace std;

void Memory_allocation_error(const char *var, const char *func);
void Cannot_open_file(const char *fname, const char *func);
void Cannot_open_file_but_continue(const char *fname, const char *func);
