/***************************************************************************

colour.h

Header file for Colour classes

Rewrite of Python Raytracer in C++

***************************************************************************/

#ifndef COLOUR_H
#define COLOUR_H

#define COLOURSPACE 3

#include <memory>

class RGB {

	public:
	
		double value[COLOURSPACE];
		
		RGB();
		RGB(const RGB &);
		RGB(double);
		RGB(double, double, double);
		
		RGB operator*(double);
		RGB operator*(const RGB);
		RGB operator+(const RGB);
		
		RGB inverse();
		std::string repr();
};

typedef std::unique_ptr<RGB> RGBPtr;

#endif
