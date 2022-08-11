// To compile, use Microsoft Visual Studio C++ compiler
//
// Header file for utils.cpp

#pragma once

void WriteMessage(char *str);
void WriteMessage(stringstream &str);
bool IsFileExist(char *fname);
void CloseProgramm(int err_code);
void StopIfErrorReturn(int err_code,char *FuncName);

