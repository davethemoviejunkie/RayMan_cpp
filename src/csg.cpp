/***************************************************************************

csg.cpp

Constructive Solid Geometry definitions etc

Rewrite of Python Raytracer in C++

***************************************************************************/

#include <iostream>

#include "csg.h"

CSGnode::CSGnode() {

	// default constructor
	
	mode = UNSET;
	left = nullptr;
	right = nullptr;
	prim = NULL;

}

CSGnode::~CSGnode() {

	std::cout << "CSGnode deleted\n";
	
}

CSGnode::CSGnode(PrimPtr p) {

	// Make a Primitive node
	
	mode = PRIMITIVE;
	left = nullptr;
	right = nullptr;
	prim = std::move(p);

}

void CSGnode::shift(double x, double y, double z) {

	// apply a shift down through all the primitives
	
	if (prim)
		prim->shift(x, y, z);
	else {
		left->shift(x, y, z);
		right->shift(x, y, z);
	}

}

void CSGnode::stretch(double s) {

	// apply a stretch down through all the primitives
	
	if (prim)
		prim->stretch(s);
	else {
		left->stretch(s);
		right->stretch(s);
	}

}

void CSGnode::stretch(double x, double y, double z) {

	// apply a shift down through all the primitives
	
	if (prim)
		prim->stretch(x, y, z);
	else {
		left->stretch(x, y, z);
		right->stretch(x, y, z);
	}

}

void CSGnode::rotate(Axis a, double theta) {

	// apply a rotation down through all the primitives
	
	if (prim)
		prim->rotate(a, theta);
	else {
		left->rotate(a, theta);
		right->rotate(a, theta);
	}

}

void CSGnode::setColour(RGB c) {
	
	// set the colour of any attached primitives
	
	if (prim)
		prim->colour = c;
	else {
		left->setColour(c);
		right->setColour(c);
	}

}

void CSGnode::setTexture(std::function<RGB(Point &)> texture) {

	// set the texture of any attached prmitives

	if (prim)
		prim->getTexture = texture;
	else {
		left->setTexture(texture);
		right->setTexture(texture);
	}

}

void CSGnode::setPhong(double phong) {

	// set the phong of any attached prmitives

	if (prim)
		prim->phong = phong;
	else {
		left->setPhong(phong);
		right->setPhong(phong);
	}

}

void CSGnode::setDiffuse(RGB d) {

	// set the diffuse reflection of any attached prmitives

	if (prim)
		prim->diffuse = std::make_unique<RGB>(d);
	else {
		left->setDiffuse(d);
		right->setDiffuse(d);
	}

}

void CSGnode::setOpacity(RGB d) {

	// set the opacity of any attached prmitives

	if (prim)
		prim->opacity = std::make_unique<RGB>(d);
	else {
		left->setOpacity(d);
		right->setOpacity(d);
	}

}

void CSGnode::setOpacity(RGB d, double ri) {

	// set the opacity and refractive index of any attached prmitives

	if (prim) {
		prim->refractiveIndex = ri;
		prim->opacity = std::make_unique<RGB>(d);
	}
	else {
		left->setOpacity(d, ri);
		right->setOpacity(d, ri);
	}

}

void showiList(IntList i) {
	
	IntList::iterator thisInt;
	for (thisInt = i.begin(); thisInt != i.end(); thisInt++)
		std::cout << "Intersect: t = " << thisInt->t << ", with " << thisInt->prim->name << "\n";
}

IntList CSGnode::getIntersections(Ray &r, bool verbose) {

	// Get intersections from Primitives left and right
	// Combine using the node type
	
	if (verbose)
		std::cout << "\nIn get Intersections ... \n";
	
	IntList iList;
	
	if (mode == PRIMITIVE) {
		// call the primitive intersection routine
		if (verbose)
			std::cout << "\nCSGnode:getIntersections: Found a primitive (" + prim->name + ")";
		iList = prim->getIntersection(r);
	}
	else {
		// get left and right intersections
		if (verbose) {
			std::cout << "\nCSGnode:getIntersections: Getting left and right intersections\n";
			if (left) {
				std::cout << "Looking at a ";
				std::cout << ((left->prim) ? "Primitive (" + left->prim->name + ")" : "Node");
				std::cout << " on the left and a ";
				std::cout << ((right->prim) ? "Primitive (" + right->prim->name + ")" : "Node");
				std::cout << " on the right\n";
			}
		}
		IntList 
			leftSide = left->getIntersections(r, verbose),
			rightSide = right->getIntersections(r, verbose);
		if (verbose) {
			std::cout << "Left side:\n";
			showiList(leftSide);
			std::cout << "Right side:\n";
			showiList(rightSide);
		}
		// so now build up list of entering/leaving solids based on CSG
		// Assume pairs of entry/exits (apart from last which may be an entry only)
		
		if (!rightSide.size()) {
			if ((mode == ADD) || (mode == SUB))
				iList = leftSide; // adding/subtracting nothing leaves just the left
			// Note: iList is already empty, intersection with nothing is nothing
		}
		else if (!leftSide.size()) {
			if (mode == ADD)
				iList = rightSide; // Adding right to zip is just right
			// Again, subtracting from nothing or intersecting with nothing is nothing
		}
		else {
			// if adding, simply add both sides
			// subtracting, remove intersection points from left that are in the solid region of right
			// * don't forget to add the points from right "flipped" as new points
			// intersecting, return intersections points that are "inside" both left and right solid regions
			
			IntList::iterator lInt, rInt;
			int rSolid = 0, lSolid = 0, lastL = 0, lastR = 0;
			Intersection thisInt;
				
			if (mode == ADD) {
				// Simply add the two sides together
				if (verbose)
					std::cout<< "\nCSGnode:getIntersections: ADDing l + r\n";
				leftSide.merge(rightSide);
				iList = leftSide;
			}
			else {
				
				// go through intersections in order, adjusting as necessary
				if (verbose)
					std::cout << "\nCSGnode:getIntersections: " << ((mode == SUB) ? "SUBbing l - r\n" : "INTERSECTing l * r\n");
				
				while (leftSide.size() && rightSide.size()) {
					// still intersections to process
					
					lInt = leftSide.begin();
					rInt = rightSide.begin();
					if ((*rInt < *lInt) || (rInt->state == INSIDE)) {
						// process right intersection
						thisInt = rightSide.front();
						rightSide.pop_front();
						lastR = rSolid;
						rSolid += (thisInt.state == LEAVE) ? -1 : 1;
						if (verbose)
							std::cout << "Process right (t " << thisInt.t << ") lSolid " << lSolid << " rSolid " << rSolid << " state " << thisInt.state << "\n";
						if (mode == SUB) {
							// If we are inside a subtracted solid, check: have we JUST entered?
							// OR have we just LEFT a subtracted solid?
							// Either way, if there is solid on left we need to set this intersection
							if (lSolid && ((!lastR && rSolid == 1) || (lastR == 1 && !rSolid))) {
								if (verbose)
									std::cout << "Mode is SUB, Last R " << lastR << " rSolid " << rSolid << " lSolid " << lSolid << "\n";
								thisInt.normal = -thisInt.normal;
								thisInt.state = (rSolid) ? LEAVE : ENTER;
								for (int i = 0; i < lSolid; i++) // create lSolid new intersections
									iList.push_back(thisInt);
							}
							// else just quietly ignore this intersection
						}
						else if (mode == INTERSECT) {
							// keep this one IF we are inside solid on both sides
							if (lSolid && rSolid)
								iList.push_back(thisInt);
						}
						else {
							std::cout << "Unknown mode in Intersect!\n";
							exit(1);
						}
					}
					else {
						// process left intersection
						thisInt = leftSide.front();
						leftSide.pop_front();
						lastL = lSolid;
						lSolid += (thisInt.state == LEAVE) ? -1 : 1;
						if (verbose)
							std::cout << "Process left (t " << thisInt.t << ") lSolid " << lSolid << " rSolid " << rSolid << " state " << thisInt.state << "\n";
						if (mode == SUB) {
							// If we are inside a subtracted solid, leave this one alone, otherwise add it in
							if (!rSolid) {
								if (verbose)
									std::cout << "Not inside subtracted solid, so include\n";
								iList.push_back(thisInt);
							}
							else {
								if (verbose)
									std::cout << "Inside subtracted solid, so ignore\n";
							}
						}
						else if (mode == INTERSECT) {
							// keep this one IF we are inside solid on both sides
							if (lSolid && rSolid)
								iList.push_back(thisInt);
						}
						else {
							std::cout << "Unknown mode in Intersect!\n";
							exit(1);
						}						
					}
				}
				if (verbose)
					std::cout << "One side has finished ...\n";
				if (leftSide.size()) {
					// all right has gone, only left remaining
					if (verbose)
						std::cout << "Only left remaining\n";
					// if subtracting, check: are we in an infinite solid? If so, there is no more left
					// otherwise, we take the remaining left and append it
					if (mode == SUB) {
						if (!rSolid)
							iList.splice(iList.end(), leftSide);
					}
					else if (mode == INTERSECT) {
						// check -- do we still have solid on right? If so, iterate through the remaining left doing what we were doing above ...
						if (rSolid) {
							for (lInt = leftSide.begin(); lInt != leftSide.end(); lInt++) {
								lSolid += (lInt->state == LEAVE) ? -1 : 1;
								if (lSolid)
									iList.push_back(thisInt);
							}
						}
					}
					else {
						std::cout << "Uknown mode in Intersect (no right)!\n";
						exit(1);
					}
				}
				else if (rightSide.size()){
					// all left has gone, only right remaining
					if (verbose)
						std::cout << "Only right remaining\n";
					// only need to check rest of right if we have left solid to infinity ...
					if (lSolid) {
						if (verbose)
							std::cout << "Solid on left to inifinity ...\n";
						for (rInt = rightSide.begin(); rInt != rightSide.end(); rInt++) {
							lastR = rSolid;
							rSolid += (rInt->state == LEAVE) ? -1 : 1;
							if (verbose)
								std::cout <<"Next right intersection gives rSolid " << rSolid << " and lastR " << lastR << "\n";
							if (mode == SUB && ((!lastR && rSolid == 1) || (lastR == 1 && !rSolid))) {
								if (verbose)
									std::cout << "Include this intersection (flipped)\n";
								thisInt = *rInt;
								thisInt.normal = -thisInt.normal;
								thisInt.state = (rSolid) ? LEAVE : ENTER;
								for (int i = 0; i < lSolid; i++)
									iList.push_back(thisInt);
							}
							else if (mode == INTERSECT && rSolid)
								iList.push_back(*rInt);
						}
					}
				}
			}
		}
	}
	
	if (verbose) {
		std::cout << "\n\nReturning from intersect:\n";
		printIntersections(iList);
	}
	return iList;

}

CSGPtr combine(CSGPtr l, CSGnodeType m, CSGPtr r) {
	
	// combine two sub-trees into a larger tree
	
	CSGPtr newRoot = std::make_unique<CSGnode>();
	newRoot->mode = m;
	newRoot->left = std::move(l);
	newRoot->right = std::move(r);
	return newRoot;
}

CSGPtr combine(PrimPtr l, CSGnodeType m, PrimPtr r) {
	
	// combine two Primitives into a CSG tree
	
	CSGPtr newRoot = std::make_unique<CSGnode>();
	newRoot->left = std::make_unique<CSGnode>(std::move(l));
	newRoot->right = std::make_unique<CSGnode>(std::move(r));
	newRoot->mode = m;
	return newRoot;
}

CSGPtr combine(CSGPtr l, CSGnodeType m, PrimPtr r) {
	
	// graft a Primitive into a CSG tree
	
	CSGPtr newRoot = std::make_unique<CSGnode>();
	newRoot->right = std::make_unique<CSGnode>(std::move(r));
	newRoot->left = std::move(l);
	newRoot->mode = m;
	return newRoot;
}

CSGPtr combine(PrimPtr l, CSGnodeType m, CSGPtr r) {
	
	// graft a Primitive into a CSG tree
	
	CSGPtr newRoot = std::make_unique<CSGnode>();
	newRoot->mode = m;
	newRoot->left = std::make_unique<CSGnode>(std::move(l));
	newRoot->right = std::move(r);
	return newRoot;
}

void printCSGTree(CSGnode *n, int level) {

	// Display a CSG tree
	
	std::string tabs = std::string(level*2, ' ');
	
	if (n) {
		
		switch(n->mode) {
			
			case PRIMITIVE:
				std::cout << tabs << "Primitive: " << n->prim->name << "\n";
				break;
			
			case ADD:
				std::cout << tabs << "Add:\n";
				printCSGTree(n->left.get(), level + 1);
				printCSGTree(n->right.get(), level + 1);
				break;

			case SUB:
				std::cout << tabs << "Sub:\n";
				printCSGTree(n->left.get(), level + 1);
				printCSGTree(n->right.get(), level + 1);
				break;

			case INTERSECT:
				std::cout << tabs << "Intersect:\n";
				printCSGTree(n->left.get(), level + 1);
				printCSGTree(n->right.get(), level + 1);
				break;
		}
	}

}
