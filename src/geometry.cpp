/***************************************************************************

geometry.cpp

Geometry Routines -- Points, Vectors and Matrices

Rewrite of Python Raytracer in C++

***************************************************************************/

#include <iostream>
#include <cmath>

#include "geometry.h"

Vector::Vector() {

	// Default initialiser -- (0, 0, 0)
	
	for (int i = 0; i < MAXD; i++)
		coord[i] = 0.0;

}

Vector::Vector(double x, double y, double z) {

	// Initialise a vector to (x, y, z)
	
	coord[0] = x;
	coord[1] = y;
	coord[2] = z;

}



std::string Vector::repr() {
	 
	// Return a string representing this vector
	
	std::string result = "Vector(";
	for (int i = 0; i < MAXD; i++) {
		result += std::to_string(coord[i]);
		if (i < (MAXD - 1))
			result += ", ";
	}
	
	return result + ")";
	 
}

Vector Vector::operator+(const Vector &v) {
	
	// Add two vectors
	
	Vector result;
	for (int i = 0; i < MAXD; i++) 
		result.coord[i] = coord[i] + v.coord[i];
	return result;
}

Vector Vector::operator-(const Vector &v) {
	
	// Subtract two vectors
	
	Vector result;
	for (int i = 0; i < MAXD; i++) 
		result.coord[i] = coord[i] - v.coord[i];
	return result;
}

Vector Vector::operator*(double s) {

	// Scale a vector by s
	
	Vector result;
	for (int i = 0; i < MAXD; i++) 
		result.coord[i] = coord[i] * s;
	return result;	

}

Vector Vector::operator*(const Vector &v) {

	// Cross-product of two vectors
	
	Vector result;
	for (int i = 0; i < MAXD; i++)
		result.coord[i] = coord[(i + 1) % MAXD] * v.coord[(i + 2) % MAXD] - coord[(i + 2) % MAXD] * v.coord[(i + 1) % MAXD];
	return result;

}

Vector Vector::operator-() {
	
	// negate vector
	Vector result;
	for (int i = 0; i < MAXD; i++)
		result.coord[i] = -coord[i];
	return result;
}

double Vector::dot(const Vector &v) {

	// Dot product of two vectors
	
	double dotProd = 0;
	for (int i = 0; i < MAXD; i++)
		dotProd += coord[i] * v.coord[i];
	return dotProd;

}

double Vector::mag() {

	// Magnitude of a vector
	
	double dotProd = 0;
	for (int i = 0; i < MAXD; i++)
		dotProd += coord[i] * coord[i];
	return std::sqrt(dotProd);

}

Vector Vector::normalised() {

	// Return a normalised vector (magnitude of 1)
	
	double m = mag();
	Vector result;
	if (m > 0)
		for (int i = 0; i < MAXD; i++)
			result.coord[i] = coord[i] / m;
	return result;

}

Point::Point() {

	// Default initialiser -- (0, 0, 0)
	
	for (int i = 0; i < MAXD; i++)
		coord[i] = 0.0;

}

Point::Point(double x, double y, double z) {

	// Initialise a point to (x, y, z)
	
	coord[0] = x;
	coord[1] = y;
	coord[2] = z;

}



std::string Point::repr() {
	 
	// Return a string representing this point
	
	std::string result = "Point(";
	for (int i = 0; i < MAXD; i++) {
		result += std::to_string(coord[i]);
		if (i < (MAXD - 1))
			result += ", ";
	}
	
	return result + ")";
	 
}

Vector Point::operator>(const Point &p2) {
	
	//Return a vector from this point to p2
	
	Vector result;
	for (int i = 0; i < MAXD; i++)
		result.coord[i] = p2.coord[i] - coord[i];
	return result;
	
}

Point Point::operator+(const Vector &v) {

	// Move a point by a vector
	
	Point result;
	for (int i = 0; i < MAXD; i++)
		result.coord[i] = coord[i] + v.coord[i];
	return result;

}

Point Point::operator-(const Vector &v) {

	// Move a point by a vector
	
	Point result;
	for (int i = 0; i < MAXD; i++)
		result.coord[i] = coord[i] - v.coord[i];
	return result;

}

Matrix::Matrix() {
	
	
	// Dummy as already initialised in Class declaration	
}

Matrix::Matrix(double values[]) {

	for (int row = 0; row <= MAXD; row ++)
		for (int col = 0; col <= MAXD; col++)
			mat[row][col] = values[row * (MAXD + 1) + col];
	
}

std::string Matrix::repr() {

	std::string result = "";
	for (int row = 0; row <= MAXD; row ++) {
		result += "[ ";
		for (int col = 0; col <= MAXD; col++)
			result += std::to_string(mat[row][col]) + " ";
		result += "]\n";
	}
	return result;

}

Matrix Matrix::operator*(double s) {
	
	// Scalar multiplication
	
	Matrix result;
	for (int row = 0; row <= MAXD; row ++)
		for (int col = 0; col <= MAXD; col++)
			result.mat[row][col] = mat[row][col] * s;
	return result;

}


Matrix Matrix::operator+(const Matrix &m) {
	
	// Add two matrices
	
	Matrix result;
	for (int row = 0; row <= MAXD; row ++)
		for (int col = 0; col <= MAXD; col++)
			result.mat[row][col] = mat[row][col] + m.mat[row][col];
	return result;

}

Matrix Matrix::operator*(const Matrix &m) {
	
	// Multiply two matrices
	
	Matrix result;
	for (int row = 0; row <= MAXD; row ++)
		for (int col = 0; col <= MAXD; col++) {
			result.mat[row][col] = 0.0;
			for (int i = 0; i <= MAXD; i++)
				result.mat[row][col] += mat[row][i] * m.mat[i][col];
		}
	return result;

}

Point Matrix::operator*(const Point &p) {

	// Transform a point by a matrix
	// Turn into a 4D coord by adding "1" as the last coord
	
	double temp[MAXD + 1];
	int i;
	for (i = 0; i < MAXD; i++)
		temp[i] = p.coord[i];
	temp[MAXD] = 1.0; 
	
	Point result;
	for (int row = 0; row < MAXD; row ++)
		for (int col = 0; col <= MAXD; col++)
			result.coord[row] += mat[row][col] * temp[col];

	return result;

}

Vector Matrix::operator*(const Vector &v) {

	// Transform a vector by a matrix
	// Turn into a 4D coord by adding "1" as the last coord
	
	double temp[MAXD + 1];
	int i;
	for (i = 0; i < MAXD; i++)
		temp[i] = v.coord[i];
	temp[MAXD] = 0.0; // Vectors do not need to be "shifted"
	
	Vector result;
	for (int row = 0; row < MAXD; row ++)
		for (int col = 0; col <= MAXD; col++)
			result.coord[row] += mat[row][col] * temp[col];

	return result;

}

Matrix Matrix::transpose() {

	// Return a transposed matrix
	
	Matrix result;
	for (int row = 0; row <= MAXD; row++)
		for (int col = 0; col <= MAXD; col++)
			result.mat[row][col] = mat[col][row];
	return result;
}

Ray::Ray() {

	// Default constructor
	
	p = Point(0, 0, 0);
	v = Vector(0, 0, 0);

}

Ray::Ray(const Point &point, const Vector &vec) {

	// Create an instance of a ray
	
	p = point;
	v = vec;

}

std::string Ray::repr() {
	
	// String representation
	std::string outline = "Ray(";
	outline += p.repr() + ", " + v.repr() + ")";
	return outline;
}


