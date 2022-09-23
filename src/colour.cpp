/***************************************************************************

colour.cpp

Colour routines -- combining colours and intensities

Rewrite of Python Raytracer in C++

***************************************************************************/

#include <iostream>

#include "colour.h"

RGB::RGB() {

	// Default constructor - black or zero intensity
	
	for (int i = 0; i < COLOURSPACE; i++)
		value[i] = 0.0;

}

RGB::RGB(const RGB &rgb) {
	
	// make a new copy on an existing RGB value
	
	for (int i = 0; i < COLOURSPACE; i++)
		value[i] = rgb.value[i];
}

RGB::RGB(double intensity) {

	// Constructor - provide a greyscale / uniform intensity
	
	for (int i = 0; i < COLOURSPACE; i++)
		value[i] = intensity;

}

RGB::RGB(double r, double g, double b) {

	// Constructor - provide a colour / intensity
	
	value[0] = r;
	value[1] = g;
	value[2] = b;

}

RGB RGB::operator*(double scale) {
	
	// scale an RGB value by a value
	
	RGB newValue;
	
	for (int i = 0; i < COLOURSPACE; i++)
		newValue.value[i] = std::max(0.0, std::min(1.0, value[i] * scale));
	
	return newValue;
}

RGB RGB::operator*(const RGB scale) {
	
	// scale an RGB value by an intensity value
	
	RGB newValue;
	
	for (int i = 0; i < COLOURSPACE; i++)
		newValue.value[i] = std::max(0.0, std::min(1.0, value[i] * scale.value[i]));
	
	return newValue;
}

RGB RGB::operator+(const RGB scale) {
	
	// add two RGB values
	
	RGB newValue;
	
	for (int i = 0; i < COLOURSPACE; i++)
		newValue.value[i] = std::max(0.0, std::min(1.0, value[i] + scale.value[i]));
	
	return newValue;
}

RGB RGB::inverse() {
	
	// provide the "opposite" RGB value
	
	RGB newValue;
	
	for (int i = 0; i < COLOURSPACE; i++)
		newValue.value[i] = 1 - value[i];
	
	return newValue;
}

std::string RGB::repr() {
	
	// provide a representation of the colour value
	
	std::string outline = "RGB(";
	for (int i = 0; i < COLOURSPACE; i++) {
		if (i)
			outline += ", ";
		outline += std::to_string(value[i]);
	}
	outline += ")";
	return outline;
}
