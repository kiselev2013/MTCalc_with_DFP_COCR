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
 *  Class for 3D-point
 *
 *  Written by Prof. Marina G. Persova
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  mpersova@mail.ru
 *  Version 1.5 December 17, 2020
*/

#pragma once

#include <math.h>

#define EPS 1e-6

namespace pv
{

class Point3D
{
private:
protected:
	double xcoord, ycoord, zcoord; // x, y and z coordinate of point
public:
//-----------------------------------------------------------
// Constructor
//-----------------------------------------------------------
	Point3D(){xcoord=ycoord=zcoord=0.0;}
//-----------------------------------------------------------
// Constructor with parameters
//-----------------------------------------------------------
	Point3D(const double &xc, const double &yc, const double &zc) 
	{ 
		xcoord=xc; ycoord=yc; zcoord=zc; 
	}
//-----------------------------------------------------------
// Destructor
//-----------------------------------------------------------
	virtual ~Point3D() {}
//-----------------------------------------------------------
// Return a reference to the x-coordinate of a point
//-----------------------------------------------------------
	double& x() { return xcoord; }
//-----------------------------------------------------------
// Return a reference to the y-coordinate of a point
//-----------------------------------------------------------
	double& y() { return ycoord; }
//-----------------------------------------------------------
// Return a reference to the z-coordinate of a point
//-----------------------------------------------------------
	double& z() { return zcoord; }
//-----------------------------------------------------------
// Return a value to the x-coordinate of a point
//-----------------------------------------------------------
	double getx() const { return xcoord; }
//-----------------------------------------------------------
// Return a value to the y-coordinate of a point
//-----------------------------------------------------------
	double gety() const { return ycoord; }
//-----------------------------------------------------------
// Return a value to the z-coordinate of a point
//-----------------------------------------------------------
	double getz() const { return zcoord; }
//-----------------------------------------------------------
// Assignment operator
//-----------------------------------------------------------
	Point3D operator=(const Point3D &p)
	{
		xcoord=p.getx();
		ycoord=p.gety();
		zcoord=p.getz();
		return *this;
	}
//-----------------------------------------------------------
// Assignment operator for scalar value
//-----------------------------------------------------------
	Point3D operator=(const double &v)
	{
		xcoord=v;
		ycoord=v;
		zcoord=v;
		return *this;
	}
//-----------------------------------------------------------
// Addition operator
//-----------------------------------------------------------
	Point3D operator+(const Point3D &p) const
	{ 
		return Point3D(xcoord+p.getx(), ycoord+p.gety(), zcoord+p.getz());
	}
//-----------------------------------------------------------
// Subtraction operator
//-----------------------------------------------------------
	Point3D operator-(const Point3D &p) const
	{ 
		return Point3D(xcoord-p.getx(), ycoord-p.gety(), zcoord-p.getz());
	}
//-----------------------------------------------------------
// Addition with assignment operator
//-----------------------------------------------------------
	Point3D operator+=(const Point3D &p)
	{ 
		xcoord+=p.getx(); ycoord+=p.gety(); zcoord+=p.getz();
		return *this;
	}
//-----------------------------------------------------------
// Subtraction with assignment operator
//-----------------------------------------------------------
	Point3D operator-=(const Point3D &p)
	{ 
		xcoord-=p.getx(); ycoord-=p.gety(); zcoord-=p.getz();
		return *this;
	}
//-----------------------------------------------------------
// Divide with assignment operator
//-----------------------------------------------------------
	Point3D operator/=(const double &v)
	{ 
		xcoord/=v; ycoord/=v; zcoord/=v;
		return *this;
	}
//-----------------------------------------------------------
// Multiplication operator for scalar value
//-----------------------------------------------------------
	Point3D operator*(const double &a) const
	{
		return Point3D(xcoord*a, ycoord*a, zcoord*a);
	}
//-----------------------------------------------------------
// Return a reference with indexing operator
//-----------------------------------------------------------
	double& operator[](const int &i)
	{
		if (i==0) return xcoord;
		if (i==1) return ycoord;
		return zcoord;
	}
//-----------------------------------------------------------
// Return a value with indexing operator
//-----------------------------------------------------------
	const double& operator[](const int &i) const
	{
		if (i==0) return xcoord;
		if (i==1) return ycoord;
		return zcoord;
	}
//-----------------------------------------------------------
// Equality operator
//-----------------------------------------------------------
	bool operator==(const Point3D &P) const
	{
		return (fabs(xcoord-P.getx())<EPS&&
			fabs(ycoord-P.gety())<EPS&&
			fabs(zcoord-P.getz())<EPS);
	}
//-----------------------------------------------------------
// Operator of text-mode reading
//-----------------------------------------------------------
	friend ofstream& operator<<(ofstream &file, const Point3D &P) 
	{
		file << P.xcoord << " " << P.ycoord << " " << P.zcoord;
		return file;
	}
//-----------------------------------------------------------
// Operator of text-mode writing
//-----------------------------------------------------------
	friend ifstream& operator>>(ifstream &file, Point3D &P)
	{
		file >> P.xcoord >> P.ycoord >> P.zcoord;
		return file;
	}
//-----------------------------------------------------------
// To determine the minimum point
//-----------------------------------------------------------
	friend Point3D minp(const Point3D &P1, const Point3D &P2)
	{
		return Point3D(min(P1.xcoord, P2.xcoord),
					   min(P1.ycoord, P2.ycoord),
					   min(P1.zcoord, P2.zcoord));
	}
//-----------------------------------------------------------
// To determine the maximum point
//-----------------------------------------------------------
	friend Point3D maxp(const Point3D &P1, const Point3D &P2)
	{
		return Point3D(max(P1.xcoord, P2.xcoord),
					   max(P1.ycoord, P2.ycoord),
					   max(P1.zcoord, P2.zcoord));
	}
//-----------------------------------------------------------
// Ñomparison less operator
//-----------------------------------------------------------
	friend bool operator<(const Point3D &P1, const Point3D &P2)
	{
		return (P1.xcoord<P2.xcoord && P1.ycoord<P2.ycoord && P1.zcoord<P2.zcoord);
	}
};

}
