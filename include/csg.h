/***************************************************************************

csg.h

Header file for CSG trees etc

Rewrite of Python Raytracer in C++

***************************************************************************/

#ifndef CSG_H
#define CSG_H

#include <functional>
#include <memory>

#include "prim.h"
#include "geometry.h"
#include "colour.h"

//#define DEBUG


enum CSGnodeType {
	PRIMITIVE,
	ADD,
	SUB,
	INTERSECT,
	UNSET
};

class CSGnode;

typedef std::unique_ptr<CSGnode> CSGPtr;

class CSGnode {
	
	public:
		CSGnodeType mode;
		//CSGnode *left, *right;
		CSGPtr left = nullptr, right = nullptr;
		PrimPtr prim;
		
		CSGnode();
		~CSGnode();
		CSGnode(PrimPtr);
		
		// these are Primitive methods, but can be applied to a node to propagate
		void shift(double, double, double);
		void stretch(double);
		void stretch(double, double, double);
		void rotate(Axis, double);
		
		// set attributes of any primitives attached to this node (or sub-nodes)
		void setColour(RGB);
		void setTexture(std::function<RGB(Point &)>);
		void setPhong(double);
		void setDiffuse(RGB);
		void setOpacity(RGB);
		void setOpacity(RGB, double);
		
		IntList getIntersections(Ray &, bool = false);

};

// Functions to create trees from either existing subtrees, or creating new ones from primitives
// Before using std::ptrs, these were CSGnode *combine's
CSGPtr combine(CSGPtr, CSGnodeType, CSGPtr);
CSGPtr combine(PrimPtr, CSGnodeType, PrimPtr);
CSGPtr combine(CSGPtr, CSGnodeType, PrimPtr);
CSGPtr combine(PrimPtr, CSGnodeType, CSGPtr);

void printCSGTree(CSGnode *, int level = 0);

#endif
