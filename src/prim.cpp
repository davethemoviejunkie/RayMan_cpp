/***************************************************************************

prim.cpp

Primitive routines -- definitions of Intersections etc

Rewrite of Python Raytracer in C++

***************************************************************************/

#include <iostream>
#include <cmath>

#include "prim.h"

//#define DEBUG_PRIM

Intersection::Intersection() {
	
	// Default constructor
	
	t = 0.0;
	point = Point(0, 0, 0);
	normal = Vector(0, 0, 0);
	//prim = NULL;
	state = UNKNOWN;
	//next = NULL;
	
}

Intersection::Intersection(double d, Point &p, Vector &v, Primitive *ob, IntersectionState s) {
	
	// ob is a Primitive * as we never CREATE a primitive here, just REFER to one
	t = d;
	point = p;
	normal = v;
	prim = ob;
	state = s;
	//next = NULL;
	
}

std::string Intersection::repr() {
	
	std::string result = "Intersection(";
	result += std::to_string(t) + ", ";
	result += point.repr() + ", ";
	result += normal.repr() + ", ";
	result += prim->name + ", ";
	switch (state) {
		case INSIDE:
			result += "INSIDE";
			break;
		case ENTER:
			result += "ENTER";
			break;
		case LEAVE:
			result += "LEAVE";
			break;
	}
	result += ")";
	return result;
}

bool Intersection::operator<(Intersection &i) {
	
	return t < i.t;
}

void printIntersections(IntList node) {
	
	IntList::iterator thisInt;
	for (thisInt = node.begin(); thisInt != node.end(); thisInt++) {
		std::cout << "Intersect(" << thisInt->t << ": " << thisInt->point.repr() << ", ";
		std::cout << thisInt->normal.repr() << ", " << thisInt->prim->name;
		std::cout << ", State: " << thisInt->state << ")\n";
	}
}

IntList planeIntersect(Primitive *prim, Ray &r) {

	// Find whether a ray intersects a plane
		
	Ray primRay = prim->toPrimitivespace(r);
	
	IntList iList;
	Intersection thisInt;
	if (primRay.p.coord[2] * primRay.v.coord[2] < 0) {
		
		// we have an intersection, either looking down or looking up
		
		float 
			d = primRay.p.coord[2] / (-primRay.v.coord[2]), // distance to intersection
			x = primRay.p.coord[0] + d * primRay.v.coord[0],
			y = primRay.p.coord[1] + d * primRay.v.coord[1];
		
		Point iPoint(x, y, 0);
		Vector n(0.0, 0.0, 1.0);
		
		
		if (primRay.v.coord[2] > 0) {
			
			// Loooking up, so below the plane to start with: need to add an Inside point
			Vector insideNormal; // normal for an "inside" point is meaningless
			thisInt = prim->toWorldspace(0.0, primRay.p, insideNormal, INSIDE);
			iList.push_back(thisInt);
			thisInt = prim->toWorldspace(d, iPoint, n, LEAVE);
			iList.push_back(thisInt);
		}
		else {
			// Looking down so only the entry point
			thisInt = prim->toWorldspace(d, iPoint, n, ENTER);
			iList.push_back(thisInt);
		}
		
	}
	else {
		
		// No intersection with plane boundary in front of us, BUT
		// we MIGHT still be in the plane. Return an "inside" if we are ...
		if (primRay.p.coord[2] < 0) {
			Vector insideNormal; // normal for an "inside" point is meaningless
			thisInt = prim->toWorldspace(0.0, primRay.p, insideNormal, INSIDE);
			iList.push_back(thisInt);
		}
	}
	
	return iList;

}

IntList sphereIntersect(Primitive *prim, Ray &r) {
	
	// Find whether a ray intersects a sphere
	
	Ray primRay = prim->toPrimitivespace(r);
	#ifdef DEBUG_PRIM
	std::cout << "WorldRay " << r.repr() << "\nPrimRay " << primRay.repr() << "\n";
	#endif
	Point hitPt;
	Vector norm;
	IntList iList;
	Intersection thisInt;
	
	double
		A = pow(primRay.v.coord[0], 2.0) + pow(primRay.v.coord[1], 2.0) + pow(primRay.v.coord[2], 2.0),
		B = 2 * (primRay.v.coord[0] * primRay.p.coord[0] + primRay.v.coord[1] * primRay.p.coord[1] + primRay.v.coord[2] * primRay.p.coord[2]),
		C = pow(primRay.p.coord[0], 2.0) + pow(primRay.p.coord[1], 2.0) + pow(primRay.p.coord[2], 2.0) - 1; // 1 is radius of sphere squared
		
	double
		m = -B/A/2.0,
		d_sq = pow(m, 2.0) - C/A,
		t, t0, t1;
	
	if (d_sq > 0) {
		t0 = m - std::sqrt(d_sq);
		t1 = m + std::sqrt(d_sq);
		t = std::min(t0, t1);
		if (t > 0) {
			
			// Two distinct intersections
			hitPt = primRay.p + (primRay.v * t);
			norm = Point(0.0, 0.0, 0.0) > hitPt;
			#ifdef DEBUG_PRIM
			std::cout << "PRIM: sphere hitPt " << hitPt.repr() << " norm " << norm.repr() << "\n";
			#endif
			thisInt = prim->toWorldspace(t, hitPt, norm, ENTER);
			iList.push_back(thisInt);
			t = std::max(t0, t1);
			hitPt = primRay.p + (primRay.v * t);
			norm = Point(0.0, 0.0, 0.0) > hitPt;
			thisInt = prim->toWorldspace(t, hitPt, norm, LEAVE);
			iList.push_back(thisInt);
		}
		else {
			t = std::max(t0, t1);
			if (t > 0) {
				// Only one intersection in front of us, so we must start inside the sphere
				hitPt = primRay.p;
				norm = Vector(0.0, 0.0, 0.0);
				#ifdef DEBUG_PRIM
				std::cout << "PRIM: sphere single intersection, inside\n";
				#endif
				thisInt = prim->toWorldspace(t, hitPt, norm, INSIDE);
				iList.push_back(thisInt);
				hitPt = primRay.p + (primRay.v * t);
				norm = Point(0.0, 0.0, 0.0) > hitPt;
				thisInt = prim->toWorldspace(t, hitPt, norm, LEAVE);
				iList.push_back(thisInt);
			}
			else {
				#ifdef DEBUG_PRIM
				std::cout << "PRIM: sphere no hits\n";
				#endif
			}
		}
	}
	else {
		#ifdef DEBUG_PRIM
		std::cout << "PRIM: sphere missed entirely\n";
		#endif
	}
	
	return iList;

}

IntList cylIntersect(Primitive *prim, Ray &r) {

	// Find whether a ray intersects a cylinder
		
	Ray primRay = prim->toPrimitivespace(r);
	Point hitPt;
	Vector norm;
	IntList iList;
	Intersection thisInt;
	
	double
		A = pow(primRay.v.coord[0], 2.0) + pow(primRay.v.coord[1], 2.0),
		B = 2 * (primRay.v.coord[0] * primRay.p.coord[0] + primRay.v.coord[1] * primRay.p.coord[1]),
		C = pow(primRay.p.coord[0], 2.0) + pow(primRay.p.coord[1], 2.0) - 1; // 1 is radius of cylinder squared
		
	double
		m = -B/A/2.0,
		d_sq = pow(m, 2.0) - C/A,
		t, t0, t1;
	
	if (d_sq > 0) {
		t0 = m - std::sqrt(d_sq);
		t1 = m + std::sqrt(d_sq);
		t = std::min(t0, t1);
		if (t > 0) {
			
			// Two distinct intersections
			hitPt = primRay.p + (primRay.v * t);
			norm = Point(0.0, 0.0, hitPt.coord[2]) > hitPt;
			thisInt = prim->toWorldspace(t, hitPt, norm, ENTER);
			iList.push_back(thisInt);
			t = std::max(t0, t1);
			hitPt = primRay.p + (primRay.v * t);
			norm = Point(0.0, 0.0, hitPt.coord[2]) > hitPt;
			thisInt = prim->toWorldspace(t, hitPt, norm, LEAVE);
			iList.push_back(thisInt);
		}
		else {
			
			t = std::max(t0, t1);
			if (t > 0) {
				
				// Only one intersection in front of us, so we must start inside the cylinder
				hitPt = primRay.p;
				norm = Vector(0.0, 0.0, 0.0);
				thisInt = prim->toWorldspace(t, hitPt, norm, INSIDE);
				iList.push_back(thisInt);
				hitPt = primRay.p + (primRay.v * t);
				norm = Point(0.0, 0.0, hitPt.coord[2]) > hitPt;
				thisInt = prim->toWorldspace(t, hitPt, norm, LEAVE);
				iList.push_back(thisInt);
			}
			
		}
	}
	
	return iList;

}

IntList coneIntersect(Primitive *prim, Ray &r) {

	// Find whether a ray intersects a cone
		
	Ray primRay = prim->toPrimitivespace(r);
	Point hitPt;
	Vector norm;
	IntList iList;
	Intersection thisInt;
	
	double
		A = pow(primRay.v.coord[0], 2.0) + pow(primRay.v.coord[1], 2.0) - pow(primRay.v.coord[2], 2.0),
		B = 2 * (primRay.v.coord[0] * primRay.p.coord[0] + primRay.v.coord[1] * primRay.p.coord[1] - primRay.v.coord[2] * primRay.p.coord[2]),
		C = pow(primRay.p.coord[0], 2.0) + pow(primRay.p.coord[1], 2.0) - pow(primRay.p.coord[2], 2.0); // z is radius of cylinder
		
	double
		m = -B/A/2.0,
		d_sq = pow(m, 2.0) - C/A,
		t, t0, t1;
	
	if (d_sq > 0) {
		t0 = m - std::sqrt(d_sq);
		t1 = m + std::sqrt(d_sq);
		t = std::min(t0, t1);
		if (t > 0) {
			
			// Two distinct intersections
			hitPt = primRay.p + (primRay.v * t);
			norm = Point(0.0, 0.0, hitPt.coord[2] * 2) > hitPt;
			thisInt = prim->toWorldspace(t, hitPt, norm, ENTER);
			iList.push_back(thisInt);
			t = std::max(t0, t1);
			hitPt = primRay.p + (primRay.v * t);
			norm = Point(0.0, 0.0, hitPt.coord[2] * 2) > hitPt;
			thisInt = prim->toWorldspace(t, hitPt, norm, LEAVE);
			iList.push_back(thisInt);
		}
		else {
			
			t = std::max(t0, t1);
			if (t > 0) {
			
				// Only one intersection in front of us, so we must start inside the cylinder
				hitPt = primRay.p;
				norm = Vector(0.0, 0.0, 0.0);
				thisInt = prim->toWorldspace(t, hitPt, norm, INSIDE);
				iList.push_back(thisInt);
				hitPt = primRay.p + (primRay.v * t);
				norm = Point(0.0, 0.0, hitPt.coord[2] * 2) > hitPt;
				thisInt = prim->toWorldspace(t, hitPt, norm, LEAVE);
				iList.push_back(thisInt);
			}
		}
	}
	
	return iList;

}

Primitive::Primitive(PrimitiveType prim, std::string theName) {

	// Constructor
	
	// common initial setups
	colour = RGB(0.5, 0.5, 0.5); //mid-grey
	//diffuse = NULL; // taken care of with switch to smart pointers
	phong = 0.0;
	//opacity = NULL; // taken care of with switch to smart pointers
	refractiveIndex = 1.0;
	getTexture = nullptr;
	
	switch (prim) {
		
		case PLANE:
			name = (theName == "") ? "Plane" : theName;
			doIntersection = planeIntersect;
			break;
		
		case SPHERE:
			name = (theName == "") ? "Sphere" : theName;
			doIntersection = sphereIntersect;
			break;
		
		case CYLINDER:
			name = (theName == "") ? "Cylinder" : theName;
			doIntersection = cylIntersect;
			break;
		
		case CONE:
			name = (theName == "") ? "Cone" : theName;
			doIntersection = coneIntersect;
			break;
	}

}

Primitive::~Primitive() {

	std::cout << "Primitive " << name << " deleted.\n";

}

IntList Primitive::getIntersection(Ray &r) {

	// Standard entry point for grabbing an iList
	return doIntersection(this, r);

}

std::string Primitive::repr() {

	std::string result;
	
	result = "Primitive: " + name + "\nTransformation:\n";
	result += transformation.repr() + "Inverse:\n" + inverseTransformation.repr();
	result += "Colour: " + colour.repr() + "\n";
	return result;

}

void Primitive::apply(Matrix &transform, Matrix &inverse) {

	// Apply a transformation matrix to this primitive
	
	transformation = transform * transformation;
	inverseTransformation = inverseTransformation * inverse;

}

void Primitive::stretch(double s) {
	
	// Scale a primitive by s
	// Scale centre is origin (0, 0, 0)
	
	Matrix toWorld, toPrim;
	
	for (int i = 0; i < MAXD; i++) {
		toWorld.mat[i][i] = s;
		toPrim.mat[i][i] = 1.0 / s; // inverse transformation
	}
	apply(toWorld, toPrim);

}

void Primitive::stretch(double sx, double sy, double sz) {
	
	// Scale a primitive by (sx, sy, sz)
	// Scale centre is origin (0, 0, 0)
	
	Matrix toWorld, toPrim;
	
	toWorld.mat[0][0] = sx;
	toPrim.mat[0][0] = 1.0 / sx; // inverse transformation
	toWorld.mat[1][1] = sy;
	toPrim.mat[1][1] = 1.0 / sy; // inverse transformation
	toWorld.mat[2][2] = sz;
	toPrim.mat[2][2] = 1.0 / sz; // inverse transformation
	
	apply(toWorld, toPrim);

}

void Primitive::shift(double sx, double sy, double sz) {

	// Move a primitive by (sx, sy, sz)
	
	Matrix toWorld, toPrim;
	
	toWorld.mat[0][3] = sx;
	toWorld.mat[1][3] = sy;
	toWorld.mat[2][3] = sz;
	
	toPrim.mat[0][3] = -sx;
	toPrim.mat[1][3] = -sy;
	toPrim.mat[2][3] = -sz;
	
	apply(toWorld, toPrim);

}

void Primitive::rotate(Axis a, double theta) {

	// Rotate a primitive theta degrees about an axis
	
	double rad = theta / 180.0 * M_PI;
	
	Matrix toWorld, toPrim;
	switch(a) {
		
		case X_AXIS:
			toWorld.mat[1][1] = std::cos(rad);
			toWorld.mat[1][2] = -std::sin(rad);
			toWorld.mat[2][1] = std::sin(rad);
			toWorld.mat[2][2] = std::cos(rad);
			break;
		
		case Y_AXIS:
			toWorld.mat[0][0] = std::cos(rad);
			toWorld.mat[0][2] = std::sin(rad);
			toWorld.mat[2][0] = -std::sin(rad);
			toWorld.mat[2][2] = std::cos(rad);
			break;

		case Z_AXIS:
			toWorld.mat[0][0] = std::cos(rad);
			toWorld.mat[0][1] = -std::sin(rad);
			toWorld.mat[1][0] = std::sin(rad);
			toWorld.mat[1][1] = std::cos(rad);
			break;
	}
	
	toPrim = toWorld.transpose();
	apply(toWorld, toPrim);

}

Intersection Primitive::toWorldspace(double t, Point &primP, Vector &primN, IntersectionState s) {
	
	// transform a hitpoint from primitive to world space
	Point worldP;
	Vector worldN;
	worldP = transformation * primP;
	worldN = transformation * primN;
	worldN = worldN.normalised();
	Intersection i(t, worldP, worldN, this, s);
	return i;
}

Point Primitive::toPrimitivespace(Point &thePoint) {
	
	// transform a point from world to primitive space
	return inverseTransformation * thePoint;
}

Ray Primitive::toPrimitivespace(Ray &theRay) {

	// transform a ray from world to primitive space

	Point thePoint = inverseTransformation * theRay.p;
	Vector theVector = inverseTransformation * theRay.v;
	return Ray(thePoint, theVector);
}


