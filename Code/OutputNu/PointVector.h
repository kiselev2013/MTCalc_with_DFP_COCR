#pragma once

#include <math.h>

#define EPS 1e-6

namespace pv
{

class Point2D
{
private:
protected:
	double xcoord, ycoord;
public:
	Point2D()
	{ 
		xcoord=ycoord=0.; 
	}
	Point2D(const double &xc, const double &yc) 
	{ 
		xcoord=xc; ycoord=yc; 
	}
	virtual ~Point2D() {}
	double& x() { return xcoord; }
	double& y() { return ycoord; }
	double getx() const { return xcoord; }
	double gety() const { return ycoord; }
	Point2D operator=(const Point2D &p)
	{
		xcoord=p.getx();
		ycoord=p.gety();
		return *this;
	}
	Point2D operator+(const Point2D &p) const
	{ 
		return Point2D(xcoord+p.getx(), ycoord+p.gety());
	}
	Point2D operator-(const Point2D &p) const
	{ 
		return Point2D(xcoord-p.getx(), ycoord-p.gety());
	}
	Point2D operator+=(const Point2D &p)
	{ 
		xcoord+=p.getx(); ycoord+=p.gety();
		return *this;
	}
	Point2D operator-=(const Point2D &p)
	{ 
		xcoord-=p.getx(); ycoord-=p.gety();
		return *this;
	}
	Point2D operator*(const double &a) const
	{
		return Point2D(xcoord*a, ycoord*a);
	}
	bool operator==(const Point2D &P) const
	{
		return (fabs(xcoord-P.getx())<EPS&&
			fabs(ycoord-P.gety())<EPS);
	}
	friend ofstream& operator<<(ofstream &file, const Point2D &P) 
	{
		file << P.xcoord << P.ycoord;
		return file;
	}
	friend ifstream& operator>>(ifstream &file, Point2D &P)
	{
		file >> P.xcoord >> P.ycoord;
		return file;
	}
	void RoundToInt()
	{
		xcoord=(int)xcoord;
		ycoord=(int)ycoord;
	}

};

class Point3D
{
private:
protected:
	double xcoord, ycoord, zcoord;
public:
	Point3D()
	{ 
		xcoord=ycoord=zcoord=0.; 
	}
	Point3D(const double &xc, const double &yc, const double &zc) 
	{ 
		xcoord=xc; ycoord=yc; zcoord=zc; 
	}
	virtual ~Point3D() {}
	double& x() { return xcoord; }
	double& y() { return ycoord; }
	double& z() { return zcoord; }
	double getx() const { return xcoord; }
	double gety() const { return ycoord; }
	double getz() const { return zcoord; }
	Point3D operator=(const Point3D &p)
	{
		xcoord=p.getx();
		ycoord=p.gety();
		zcoord=p.getz();
		return *this;
	}
	Point3D operator=(const double &v)
	{
		xcoord=v;
		ycoord=v;
		zcoord=v;
		return *this;
	}
	Point3D operator=(const Point2D &p)
	{
		xcoord=p.getx();
		ycoord=p.gety();
		zcoord=0;
		return *this;
	}
	Point3D operator+(const Point3D &p) const
	{ 
		return Point3D(xcoord+p.getx(), ycoord+p.gety(), zcoord+p.getz());
	}
	Point3D operator-(const Point3D &p) const
	{ 
		return Point3D(xcoord-p.getx(), ycoord-p.gety(), zcoord-p.getz());
	}
	Point3D operator+=(const Point3D &p)
	{ 
		xcoord+=p.getx(); ycoord+=p.gety(); zcoord+=p.getz();
		return *this;
	}
	Point3D operator-=(const Point3D &p)
	{ 
		xcoord-=p.getx(); ycoord-=p.gety(); zcoord-=p.getz();
		return *this;
	}
	Point3D operator/=(const double &v)
	{ 
		xcoord/=v; ycoord/=v; zcoord/=v;
		return *this;
	}
	Point3D operator*(const double &a) const
	{
		return Point3D(xcoord*a, ycoord*a, zcoord*a);
	}
	double& operator[](const int &i)
	{
		if (i==0) return xcoord;
		if (i==1) return ycoord;
		return zcoord;
	}
	const double& operator[](const int &i) const
	{
		if (i==0) return xcoord;
		if (i==1) return ycoord;
		return zcoord;
	}
	bool operator==(const Point3D &P) const
	{
		return (fabs(xcoord-P.getx())<EPS&&
			fabs(ycoord-P.gety())<EPS&&
			fabs(zcoord-P.getz())<EPS);
	}
	friend ofstream& operator<<(ofstream &file, const Point3D &P) 
	{
		file << P.xcoord << " " << P.ycoord << " " << P.zcoord;
		return file;
	}
	friend ifstream& operator>>(ifstream &file, Point3D &P)
	{
		file >> P.xcoord >> P.ycoord >> P.zcoord;
		return file;
	}

	friend Point3D minp(const Point3D &P1, const Point3D &P2)
	{
		return Point3D(min(P1.xcoord, P2.xcoord),
					   min(P1.ycoord, P2.ycoord),
					   min(P1.zcoord, P2.zcoord));
	}
	friend Point3D maxp(const Point3D &P1, const Point3D &P2)
	{
		return Point3D(max(P1.xcoord, P2.xcoord),
					   max(P1.ycoord, P2.ycoord),
					   max(P1.zcoord, P2.zcoord));
	}

	void GetCoords(float *m)
	{
		m[0]=static_cast<float>(xcoord);
		m[1]=static_cast<float>(ycoord);
		m[2]=static_cast<float>(zcoord);
	}

	friend bool operator<(const Point3D &P1, const Point3D &P2)
	{
		return (P1.xcoord<P2.xcoord && P1.ycoord<P2.ycoord && P1.zcoord<P2.zcoord);
	}
};

class Vector: public Point3D
{
private:
protected:
public:
	Vector(): Point3D() {}
	Vector(const double &xc, const double &yc, const double &zc): Point3D(xc, yc, zc) {} 
	Vector(const Point3D &P) { xcoord=P.getx(); ycoord=P.gety(); zcoord=P.getz(); } 
	Vector(const Point2D &P) { xcoord=P.getx(); ycoord=P.gety(); zcoord=0; } 
	virtual ~Vector() {}
	Vector operator=(const Point3D &p)
	{
		dynamic_cast<Point3D*>(this)->operator=(p);
		return *this;
	}
	void operator/=(const double &d)
	{
		xcoord/=d; ycoord/=d; zcoord/=d; 
	}
	void operator*=(const double &d)
	{
		xcoord*=d; ycoord*=d; zcoord*=d; 
	}
	void normalize()
	{
		double s=this->norm();
		if (s>0.)
		{
			xcoord/=s; ycoord/=s; zcoord/=s;
		}
	}
	double norm()
	{
		return sqrt(xcoord*xcoord+ycoord*ycoord+zcoord*zcoord);
	}
	double sp(const Vector &v) const
	{
		return (xcoord*v.getx()+ycoord*v.gety()+zcoord*v.getz());
	}
	Vector vp(const Vector &v) const
	{ 
		return Vector(
			ycoord*v.getz()-zcoord*v.gety(),
			zcoord*v.getx()-xcoord*v.getz(),
			xcoord*v.gety()-ycoord*v.getx()); 
	}
};

}