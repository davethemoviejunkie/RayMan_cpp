/***************************************************************************

prim.h

Header file for Primitives

Rewrite of Python Raytracer in C++

***************************************************************************/

#ifndef PRIM_H
#define PRIM_H

#include <functional>
#include <vector>
#include <list>
#include <memory>

#include "geometry.h"
#include "colour.h"

enum PrimitiveType {
	PLANE,
	SPHERE,
	CYLINDER,
	CONE
};

// Intersection States

enum IntersectionState {
	ENTER,
	LEAVE,
	INSIDE,
	UNKNOWN
};

#define RI_WATER 1.33
#define RI_GLASS 1.52
#define RI_DIAMOND 2.417

class Intersection; 

class Primitive;

typedef std::list<Intersection> IntList;
typedef std::unique_ptr<Primitive> PrimPtr;

class Primitive {

	public:
		std::string name;
		Matrix transformation, inverseTransformation;
		IntList getIntersection(Ray &);
		//static std::function<IntList(Primitive *, Ray &)> doIntersection;
		IntList (*doIntersection)(Primitive *, Ray &);
		
		RGB 
			colour;
		
		RGBPtr 
			diffuse = nullptr,
			opacity = nullptr;
		
		double 
			refractiveIndex,
			phong;
		
		//RGB (*getTexture)(Point &);
		std::function<RGB(Point &)> getTexture;
		
		std::string repr();
		
		Primitive(PrimitiveType, std::string theName = ""); // constructor
		~Primitive();
		void apply(Matrix &, Matrix &);
		void stretch(double);
		void stretch(double, double, double);
		void shift(double, double, double);
		void rotate(Axis, double);
		
		Intersection toWorldspace(double, Point &, Vector &, IntersectionState);
		Point toPrimitivespace(Point &);
		Ray toPrimitivespace(Ray &);
		
};

class Intersection {
	
	public:
		double t;
		Point point;
		Vector normal;
		Primitive *prim = NULL;
		IntersectionState state;
		//IntPtr next = nullptr;
		
		Intersection();
		Intersection(double, Point &, Vector &, Primitive *, IntersectionState);
		
		std::string repr();
		
		bool operator<(Intersection &);
		
};

void printIntersections(IntList);

IntList
	planeIntersect(PrimPtr, Ray &),
	sphereIntersect(PrimPtr, Ray &),
	cylIntersect(PrimPtr, Ray &),
	coneIntersect(PrimPtr, Ray &);

#include "csg.h"

#endif
