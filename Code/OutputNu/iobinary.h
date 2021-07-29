#pragma once

ifstream& operator>(ifstream& file, short &id);
ofstream& operator<(ofstream& file, const short &id);

ifstream& operator>(ifstream& file, int &id);
ofstream& operator<(ofstream& file, const int &id);

ifstream& operator>(ifstream& file, long &id);
ofstream& operator<(ofstream& file, const long &id);

ifstream& operator>(ifstream& file, double &id);
ofstream& operator<(ofstream& file, const double &id);

ifstream& operator>(ifstream& file, float &id);
ofstream& operator<(ofstream& file, const float &id);

