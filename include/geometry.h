/***************************************************************************

geometry.h

Header file for Geometry -- Vectors, Points, Matrices and Rays/Intersections

Rewrite of Python Raytracer in C++

***************************************************************************/

#ifndef GEOMETRY_H
#define GEOMETRY_H

// Maximum dimensions
#define MAXD 3

// Axes for Rotation

enum Axis {
	X_AXIS,
	Y_AXIS,
	Z_AXIS
};

#include <string>

// Vector Class -- stores a 3D (x, y, z) value

class Vector {

	public:
		double coord[MAXD];
		
		Vector();
		Vector(double, double, double);
		std::string repr();
		
		Vector operator+(const Vector &);
		Vector operator-(const Vector &);
		Vector operator*(const Vector &);
		Vector operator*(double s);
		Vector operator-();
		
		double dot(const Vector &);
		double mag();
		Vector normalised();
		
};


// Point Class -- stores a 3D (x, y, z) value

class Point {

	public:
		double coord[MAXD];
		
		Point();
		Point(double, double, double);
		std::string repr();
		
		Vector operator>(const Point &);
		Point operator+(const Vector &);
		Point operator-(const Vector &);
		
};

class Matrix {

	public:
		double mat [MAXD+1][MAXD+1] = {
				{1.0, 0.0, 0.0, 0.0},
				{0.0, 1.0, 0.0, 0.0},
				{0.0, 0.0, 1.0, 0.0},
				{0.0, 0.0, 0.0, 1.0}
			};
		
		Matrix();
		Matrix(double []); // mainly for testing
		
		std::string repr();
		
		Matrix operator*(double); //Scalar multiplication
		Matrix operator+(const Matrix &); //Add matricies
		Matrix operator*(const Matrix &); // Multiply matrices
		
		Point operator*(const Point &); // Transform a point by a matrix
		Vector operator*(const Vector &); // Transform a vector by a matrix
		
		Matrix transpose();
};

class Ray {

	public:
		Point p;
		Vector v;
		
		Ray();
		Ray(const Point &, const Vector &);
		std::string repr();
		
};



#endif
