/**
 * GENERAL REMARKS
 *
 *  This code is freely available under the following conditions:
 *
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer.
 *
 *              MTCalc_with_DFP_COCR
 *  Header file for iobinary.cpp
 *
 *  Written by Prof. Marina G. Persova
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  mpersova@mail.ru
 *  Version 1.5 December 17, 2020
*/

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
