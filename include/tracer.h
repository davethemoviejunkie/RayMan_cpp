/***************************************************************************

tracer.h

Header file for Renderer

Rewrite of Python Raytracer in C++

***************************************************************************/

#ifndef TRACE_H
#define TRACE_H

#include <string>
#include <list>

#include "csg.h"
#include "colour.h"
#include "file.h"
#include "geometry.h"


#define DEBUG_TRACE

enum renderMode {
	
	RENDER,
	LIGHT,
	LIGHT_TRANS

};

class Screen {
	
	// Viewing screen
	
	public:
		int width, height;
		double rotation;
		
		Screen();
		Screen(int, int, double rot = 0.0);
		std::string repr();

};

class Light_directional {
	
	// Directional light sources (moons)
	
	public:
		RGB intensity;
		Vector direction;
		
		Light_directional(RGB &, Vector &);
		std::string repr();

};

class Light_point {

	// Point light sources (suns)
	
	public:
		RGB intensity;
		Point point;
		
		Light_point(RGB &, Point &);
		std::string repr();

};

class Scene {
	
	// scene definition
	
	public:
	
		Screen screen;
		RGB background;
		Point viewPoint, toPoint;
		CSGPtr objectTree;
		int transmissionCount;
		RGB ambient;
		std::list<Light_directional> moon;
		std::list<Light_point> sun;
		double refractiveIndex;
		int antialias;
		
		Scene();
		~Scene();
		std::string repr();
		
		void addMoon(RGB, Vector);
		void addSun(RGB, Point);
};

void render(Scene &, std::string, bool = false);
RGB trace(renderMode, Ray, Scene &, int = -1, double = -1.0, bool = false);

#endif
