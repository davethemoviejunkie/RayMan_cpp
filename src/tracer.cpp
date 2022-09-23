/***************************************************************************

tracer.cpp

Tracing routines -- casting rays and building the picture

Rewrite of Python Raytracer in C++

***************************************************************************/

#include <iostream>
#include <vector>
#include <cmath>

#include "tracer.h"


Screen::Screen (int w, int h, double r) {
	
	// default viewing screen
	
	width = w;
	height = h;
	rotation = r;
	
}

Screen::Screen () {
	
	// default constructor
	
	width = 0;
	height = 0;
	rotation = 0.0;
	
}

std::string Screen::repr() {
	
	std::string result = "Screen(";
	result += std::to_string(width) + ", " + std::to_string(height) + ", " + std::to_string(rotation) + ")";
	return result;
}

Light_directional::Light_directional(RGB &c, Vector &v) {
	
	// Object definition for a moon
	
	intensity = c;
	direction = v;
}

std::string Light_directional::repr() {

	std::string result = "Moon(";
	result += intensity.repr() + ", " + direction.repr() + ")";
	return result;
}

Light_point::Light_point(RGB &c, Point &p) {
	
	// Object definition for a sun
	
	intensity = c;
	point = p;	
}

std::string Light_point::repr() {

	std::string result = "Sun(";
	result += intensity.repr() + ", " + point.repr() + ")";
	return result;
}

Scene::Scene() {
	
	// Create a scene
	
	screen = Screen(800, 600);
	background = RGB(0.1, 0.1, 0.1);
	viewPoint = Point(500, 500, 500);
	toPoint = Point(0, 0, 0);
	objectTree = NULL;
	transmissionCount = 0;
	ambient = RGB(1, 1, 1);
	//moon = NULL;
	//sun = NULL;
	refractiveIndex = 1.0;
	antialias = 0;
	
}

Scene::~Scene() {
	
	std::cout << "Scene deleted\n";

}

std::string Scene::repr() {

	std::string result = "Scene:\n";
	result += "\t" + screen.repr() + "\n";
	result += "\tBg: " + background.repr() + "\n";
	result += "\tView: " + viewPoint.repr() + "\n";
	result += "\tTo: " + toPoint.repr() + "\n";
	result += "\tTransmission count: " + std::to_string(transmissionCount) + "\n";
	result += "\tAmbient: " + ambient.repr() + "\n";
	result += "\tLighting:\n";
	std::list<Light_directional>::iterator m;
	for (m = moon.begin(); m != moon.end(); m++)
		result += "\t\tDirectional: " + m->intensity.repr() + " " + m->direction.repr() + "\n";
	std::list<Light_point>::iterator s;
	for (s = sun.begin(); s != sun.end(); s++)
		result += "\t\tPoint: " + s->intensity.repr() + " " + s->point.repr() + "\n";
	result += "Refractive index: " + std::to_string(refractiveIndex) + "\n";
	result += "Antialias: " + std::to_string(antialias) + "\n";

	return result;
}

void Scene::addMoon(RGB intensity, Vector dir) {

	Light_directional l(intensity, dir);
	this->moon.push_back(l);

}

void Scene::addSun(RGB intensity, Point p) {

	Light_point l(intensity, p);
	this->sun.push_back(l);

}

/****************** Main Rendering routines *****************/

RGB traceLight(Ray r, Scene &s, int transmissionCount, double currentRI, std::list<Light_directional>::iterator m, bool verbose) {
	
	// tracing shadow rays
	
	if (transmissionCount == 0) {// reached the end of our recursion limit, so just return black
		if (verbose)
			std::cout << "End of recursion, returning black\n";
		return RGB(0, 0, 0);
	}

	if (transmissionCount < 0)
		transmissionCount = s.transmissionCount;
	
	if (currentRI < 0.0)
		currentRI = s.refractiveIndex;
	
	Intersection hitPt;
	
	IntList iList = s.objectTree->getIntersections(r, false);

	if (iList.size()) {
		// we've hit something! Is it opaque?
		hitPt = iList.front();
		if (!hitPt.prim->opacity) // return black
			return RGB(0, 0, 0);
			
		// Otherwise, simply turn this ray into a colour and multiply that by the moon intensity (fudge it)
		return *hitPt.prim->opacity;
		RGB tracedColour = trace(RENDER, r, s, transmissionCount - 1, currentRI, verbose);
		
		return tracedColour;
		
		// O ... K ... let's refract through and see what we can find
		
		// repeat of Snell's law from below ... can we make this more efficient?
		double
			n1 = currentRI,
			n2 = hitPt.prim->refractiveIndex;
		Vector
			rNormal = hitPt.normal;
		if (hitPt.state == LEAVE) { // we're leaving, so find the next RI and flip the normal
			if (iList.size() > 1) {
				IntList::iterator i = iList.begin();
				i++;
				if (i->state == LEAVE) { // we're leaving the next solid, so the next solid is n2
					n2 = i->prim->refractiveIndex;
				}
				else { // we're entering the next solid, so this must be empty space
					n2 = s.refractiveIndex;
				}
			}
			else { // no other objects
				n2 = s.refractiveIndex;
			}
			rNormal = -rNormal;
		}
		double
			refractiveRatio = n1 / n2,
			c1 = rNormal.dot(r.v.normalised()),
			k = 1 - (refractiveRatio * refractiveRatio) * (1 - c1 * c1);
		if (k < 0) {
			
			// we have total internal reflection
			
			// just return this object's colour
			
			return iList.front().prim->colour;
		}
		else {
			
			// let's refract
			
			Vector transmitted = (r.v * refractiveRatio) + (rNormal * (refractiveRatio * c1 - std::sqrt(k)));
			Ray transmittedRay;
			transmittedRay.p = iList.front().point - (rNormal * 0.5);
			transmittedRay.v = transmitted;
			RGB refractedColour = traceLight(transmittedRay, s, transmissionCount - 1, n2, m, verbose);
			if (iList.front().state == LEAVE) // colour this returned light by the opacity of this object 
				refractedColour = refractedColour * (*iList.front().prim->opacity);
			return refractedColour;
		}
	}
	else {
		return RGB(1, 1, 1);
		// OK, we've hit empty space
		RGB light(1, 1, 1);
		double lightScale = r.v.dot(-m->direction);
		if (verbose)
			std::cout << "Dotting " << r.v.repr() << " and " << m->direction.repr() << " gives " << lightScale << "\n";
		RGB result = light * lightScale;
		if (verbose)
			std::cout << "Found empty space, returning " << result.repr() << "\n";
		return result;
	}
}


RGB trace(renderMode mode, Ray r, Scene &s, int transmissionCount, double currentRI, bool verbose) {

	// cast a ray into a scene and return the colour value 
	
	// dummy to see if we can create an image
	
	if (transmissionCount == 0)
		return RGB(0, 0, 0);
	
	if (verbose) {
		std::cout << "In trace! Mode is ";
		switch (mode) {
			case RENDER:
				std::cout << "Render\n";
				break;
			case LIGHT:
				std::cout << "Light-seeking\n";
				break;
			case LIGHT_TRANS:
				std::cout << "Processing a light ray\n";
				break;
		}
	}
	
	if (transmissionCount < 0)
		transmissionCount = s.transmissionCount;
	
	if (currentRI < 0.0)
		currentRI = s.refractiveIndex;
	
	if (verbose)
		std::cout << "Got transmission " << transmissionCount << " and RI " << currentRI << "\n";
		
	Intersection hitPt;
	
	IntList iList = s.objectTree->getIntersections(r, verbose);
	
	if (verbose)
		if (iList.size())
			std::cout << "iList returned intersections\n";
		else
			std::cout << "iList returned no intersections\n";

	RGB shadow(0, 0, 0), directionalComponent(0, 0, 0), pixelColour(0, 0, 0), baseColour;
	Ray shadowCast;
	
	if (!iList.size()) {
		if (verbose)
			std::cout << "Hit background, returning BG colour\n";
		return s.background;
	}
	else {
		// OK, got an intersection ...
		hitPt = iList.front();
		iList.pop_front();
		// check: are we "inside" an object
		if (hitPt.state == INSIDE) {
			if (hitPt.prim->opacity) {
				// we are inside a transparent object
				// grab next intersection (if it exists)
				if (iList.size()) {
					hitPt = iList.front();
					iList.pop_front();
				}
				else
					return hitPt.prim->colour * s.ambient;
			}
			else
				return hitPt.prim->colour * s.ambient;
		}
		if (verbose) {
			std::cout << "Processing intersection " << hitPt.repr() << "\n";
		}
		
		// Not inside anything, so ...
		
		// first, grab base colour
		
		if (hitPt.prim->getTexture) {
			Point textP = hitPt.prim->transformation * hitPt.point;
			baseColour = hitPt.prim->getTexture(textP);
			//baseColour = hitPt->prim->getTexture(hitPt->point);
		}
		else {
			baseColour = hitPt.prim->colour;
		}
		if (verbose)
			std::cout << "Computing colour ... starting with " << baseColour.repr() << "\n";
		// start with ambient lighting
		
		RGB incidentLight = s.ambient;
		
		if (verbose)
			std::cout << "Light intensity ... ambient is " << incidentLight.repr() << "\n";
		
		// Now add in lights from moons ... as long as we have a normal
		// we may /not/ have a normal if we've hit the top of a cone (for instance)
		
		if (verbose)
			std::cout << "Lighting: testing a mag of " << hitPt.normal.mag() << "\n";
		
		if (hitPt.normal.mag() > 0) {
			std::list<Light_directional>::iterator m;
			shadowCast.p = hitPt.point + hitPt.normal * 0.01;
			for (m = s.moon.begin(); m != s.moon.end(); m++) {
				shadowCast.v = -(m->direction);
				if (verbose)
					std::cout << "About to cast shadow ray ... \n";
				shadow = traceLight(shadowCast, s, transmissionCount - 1, currentRI, m, verbose);
				if (verbose)
					std::cout << "Shadow returned " << shadow.repr() << "\n";
				directionalComponent = (shadow * m->intensity) * hitPt.normal.dot(shadowCast.v);
				if (verbose)
					std::cout << "Adding " << directionalComponent.repr() << " from a moon\n";
				incidentLight = incidentLight + directionalComponent;
			}
		}
		
		pixelColour = baseColour * incidentLight;
		if (verbose) {
			std::cout << "Incident light is now " << incidentLight.repr();
			std::cout << "\ngiving a pixel value of " << pixelColour.repr() << "\n";
		}
		
		// Now check for reflections/transparency
		
		if (transmissionCount && (hitPt.normal.mag() > 0)) {
		
			Vector v = r.v.normalised();
			
			if (hitPt.prim->diffuse) {
				
				// we have a reflective surface and a normal
				
				Vector reflected = v - (hitPt.normal * (2 * v.dot(hitPt.normal)));
				Ray reflectedRay;
				reflectedRay.p = hitPt.point + (reflected * 0.5);
				reflectedRay.v = reflected;
				RGB reflectedColour = trace(RENDER, reflectedRay, s, transmissionCount - 1, currentRI, verbose);
				RGB diffuse = *hitPt.prim->diffuse;
				pixelColour = (reflectedColour * diffuse) + (pixelColour * diffuse.inverse());
				
			}
			
			// now we have the "base colour", lastly check for transparency
			
			if (hitPt.prim->opacity) {
				
				// Use Snell's law to compute refracted ray, and see what we can find ...
				double n1 = currentRI, n2 = hitPt.prim->refractiveIndex; // default to entering this primitive
				Vector rNormal = hitPt.normal;
				
				// get n1, n2 and the normal sorted
				
				if (hitPt.state == LEAVE) {
					// OK, we're LEAVING the solid, so that means we need to find the NEXT solid to grab the RI of (if it exists)
					n1 = hitPt.prim->refractiveIndex; 
					if (iList.size())
						if (iList.front().state == LEAVE) // the next solid is enclosing this one
							n2 = iList.front().prim->refractiveIndex;
						else
							// We're entering a new solid and leaving this one ... let's ASSUME that there's empty space at this stage (need to add solid counts to intersections to fix this)
							n2 = s.refractiveIndex; 
					else
						// no other intersection, so we're in empty space
						n2 = s.refractiveIndex;
					rNormal = -(hitPt.normal);
				}
				
				double
					refractiveRatio = n1 / n2,
					c1 = rNormal.dot(r.v.normalised()),
					k = 1 - (refractiveRatio * refractiveRatio) * (1 - c1 * c1);
				if (k < 0) {
					
					// we have total internal reflection
					
					Vector reflected = v + (rNormal * (2 * v.dot(rNormal)));
					Ray reflectedRay;
					reflectedRay.p = hitPt.point + (reflected * 0.5);
					reflectedRay.v = reflected;
					RGB reflectedColour = trace(RENDER, reflectedRay, s, transmissionCount - 1, currentRI, verbose);
					RGB opac = *hitPt.prim->opacity;
					pixelColour = (reflectedColour * opac) + (pixelColour * opac.inverse());
				}
				else {
					
					Vector transmitted = (r.v * refractiveRatio) + (rNormal * (refractiveRatio * c1 - std::sqrt(k)));
					Ray transmittedRay;
					transmittedRay.p = hitPt.point - (rNormal * 0.5);
					transmittedRay.v = transmitted;
					RGB refractedColour = trace(RENDER, transmittedRay, s, transmissionCount - 1, n2, verbose);
					RGB opac = *hitPt.prim->opacity;
					pixelColour = (refractedColour * opac) + (pixelColour * opac.inverse());
				}
				
			}
		}
		
		// lastly, add Blinn-Phong highlights
		
		if ((hitPt.prim->phong > 0.0) && (hitPt.normal.mag() > 0.0)) {
			Vector cameraDir = -r.v;
			if (verbose) {
				std::cout << "Doing Phong with " << hitPt.prim->phong  << ". View ray is " << r.v.repr() << " and camera dir is " << cameraDir.repr() << "\n";
			}
			std::list<Light_directional>::iterator m;
			for (m = s.moon.begin(); m != s.moon.end(); m++) {
				Vector halfDirection = cameraDir - m->direction;
				halfDirection = halfDirection.normalised();
				if (verbose) {
					std::cout << "Moon dir: " << m->direction.repr() << "\nhalfDir: " << halfDirection.repr() << "\n";
				}
				double blinnPhong;
				blinnPhong = pow(std::max(halfDirection.dot(hitPt.normal), 0.0), hitPt.prim->phong);
				RGB hilite = m->intensity;
				hilite = hilite * blinnPhong;
				if (verbose)
					std::cout << "Moon colour: " << m->intensity.repr() << " and scale factor " << blinnPhong << " giving a hilite of " << hilite.repr() << "\n";
				pixelColour = pixelColour + hilite;
				if (verbose)
					std::cout << "Pixel colour now " << pixelColour.repr() << "\n";
			}
		}
		
		return pixelColour;
	}
		
}

RGB trace(renderMode mode, Ray r, Scene &s, bool verbose = false) {

	return trace(mode, r, s, -1, -1.0, verbose);
}

std::string statusDots(int n) {
	
	std::string result = "";
	for (int i = 0; i < n; i++)
		result += ".";
	return result;
}

void render(Scene &theScene, std::string fname, bool verbose) {

	// Take a scene and render it
	
	if (verbose)
		std::cout << "Render: Start\n";
	
	// moderate lighting so that the maximum is 1 at any point in the scene
	/*
	double totalLight[COLOURSPACE];
	
	for (int i = 0; i < COLOURSPACE; i++)
		totalLight[i] = theScene.ambient.value[i];
	
	Light_directional *m = theScene.moon;
	while (m) {
		for (int i = 0; i < COLOURSPACE; i++)
			totalLight[i] += m->intensity.value[i];
		m = m->next;
	}
	
	double maxLight = 0.0;
	for (int i = 0; i < COLOURSPACE; i++)
		if (totalLight[i] > maxLight)
			maxLight = totalLight[i];
	
	if (maxLight > 1.0) {
		
		// scale all the lights down so that we don't overblow the scene...
		
		for (int i = 0; i < COLOURSPACE; i++)
			theScene.ambient.value[i] = theScene.ambient.value[i] / maxLight;
			
		m = theScene.moon;
			while (m) {
				for (int i = 0; i < COLOURSPACE; i++)
					m->intensity.value[i] = m->intensity.value[i] / maxLight;
				m = m->next;
			}
		
		std::cout << "Lights scaled down, max was " << maxLight << "\n";
	}
	*/
	
	int
		w = theScene.screen.width,
		h = theScene.screen.height,
		s = w * h;

	Vector 
		z(0, 0, 1),
		viewDir = (theScene.viewPoint > theScene.toPoint).normalised(),
		horizontal = (viewDir * z).normalised(),
		vertical = (viewDir * horizontal).normalised();

	if (verbose)
		std::cout << "Render: vectors initialised\n";

	Point
		screenCentre = theScene.viewPoint + (viewDir * w),
		topLeft = (screenCentre - (horizontal * (w / 2))) - (vertical * (h / 2)),
		thisPixel;

	if (verbose)
		std::cout << "Render: Points initialised\n";
	
	std::vector<RGB> img(s);

	if (verbose)
		std::cout << "Render: Image array initialised\n";

	Ray
		theRay;
		
	if (verbose)
		std::cout << "Render: variables initialised\n";
	
	theRay.p = theScene.viewPoint;
	
	for (int down = 0; down < h; down++) {
		if (!(down % 100))
			std::cout << "\nRendering row " << down << "\n";
			for (int across = 0; across < w; across++) {
				std::cout << "\r" + statusDots((down % 100) / 2);
				thisPixel = topLeft + horizontal * across + vertical * down;
				theRay.v = (theScene.viewPoint > thisPixel).normalised();
				img[down * w + across] = trace(RENDER, theRay, theScene, verbose);
			}
	}
	
	std::cout << "\nRender complete\n";
	
	outputPicture(img, w, h, fname);
	
	std::cout << "\nImage file dumped to " + fname + "\n";

}

