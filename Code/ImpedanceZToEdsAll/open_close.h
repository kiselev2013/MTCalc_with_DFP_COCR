#pragma once

int OpenInputFile(ifstream &inf,char *fname,ios_base::open_mode mode=ios_base::in);
int OpenOutputFile(ofstream &ofp,char *fname,ios_base::open_mode mode=ios_base::out);
int CloseInputFile(ifstream &inf,char *fname);
int CloseOutputFile(ofstream &ofp,char *fname);
