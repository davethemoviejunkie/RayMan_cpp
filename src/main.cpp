// main.cpp
// main testing file

#include <iostream>
#include <cmath>

#include "geometry.h"
#include "prim.h"
#include "csg.h"
#include "colour.h"
#include "tracer.h"

#define nFrames 1

RGB chessboard(Point &p) {
	
	const int squareSize = 50;
	RGB result;
	int
		xStep = ((int)p.coord[0] / squareSize) % 2,
		yStep = ((int)p.coord[1] / squareSize) % 2;
	if ((xStep + yStep) % 2)
		result = RGB(0.2, 0.8, 0.2);
	else
		result = RGB(0.8, 0.8, 0.8);
	return result;
}

int main(int argc, char * argv[]) {

	
	std::cout << "Start of program\n";
	
	Scene s;
	std::cout << "Scene created\n";
	
	PrimPtr
		basePlane = std::make_unique<Primitive>(PLANE, "Ground"),
		cyl = std::make_unique<Primitive>(CYLINDER, "Cyl"),
		cylT = std::make_unique<Primitive>(PLANE, "CylT"),
		cylB = std::make_unique<Primitive>(PLANE, "CylB"),
		cone = std::make_unique<Primitive>(CONE, "Cone"),
		coneT = std::make_unique<Primitive>(PLANE, "ConeT"),
		coneB = std::make_unique<Primitive>(PLANE, "ConeB"),
		sp = std::make_unique<Primitive>(SPHERE, "Sphere");
	
	PrimPtr
		cubeB = std::make_unique<Primitive>(PLANE, "CubeB"),
		cubeT = std::make_unique<Primitive>(PLANE, "CubeT"),
		cubeF = std::make_unique<Primitive>(PLANE, "CubeF"),
		cubeK = std::make_unique<Primitive>(PLANE, "CubeK"),
		cubeL = std::make_unique<Primitive>(PLANE, "CubeL"),
		cubeR = std::make_unique<Primitive>(PLANE, "CubeR");
	std::cout << "Primitives created\n";
	
	//Ground
	basePlane->colour = RGB(0, 1, 0);
	basePlane->getTexture = chessboard;
	basePlane->diffuse = std::make_unique<RGB>(0.4);
	
	std::cout << "Baseplane colour and texture set\n";
	
	cubeT->shift(0, 0, 600);
	cubeF->rotate(X_AXIS, -90);
	cubeL->rotate(Y_AXIS, 90);
	cubeK->rotate(X_AXIS, 90);
	cubeK->shift(0, 600, 0);
	cubeR->rotate(Y_AXIS, -90);
	cubeR->shift(600, 0, 0);
	CSGPtr cube;
	cube = combine(std::move(cubeT), SUB, std::move(cubeB));
	cube = combine(std::move(cube), SUB, std::move(cubeL));
	cube = combine(std::move(cube), SUB, std::move(cubeR));
	cube = combine(std::move(cube), SUB, std::move(cubeF));
	cube = combine(std::move(cube), SUB, std::move(cubeK));
	cube->setColour(RGB(0, 1, 0));
	cube->setDiffuse(RGB(0.4));
	cube->setOpacity(RGB(0.9), 1.5);
	cube->setPhong(50);
	cube->shift(-1100, 1700, 100); // x -1500
	
	sp->colour = RGB(1, 0, 0);
	sp->diffuse = std::make_unique<RGB>(0.4);
	sp->stretch(300);
	sp->shift(0, 2000, 400);
	sp->phong = 100;
	
	std::cout << "Red sphere done\n";
	
	cyl->stretch(300);
	cylB->shift(0, 0, -300);
	cylT->rotate(X_AXIS, 180);
	cylT->shift(0, 0, 300);
	
	CSGPtr cylinder;
	cylinder = combine(std::move(cyl), SUB, std::move(cylB));
	cylinder = combine(std::move(cylinder), SUB, std::move(cylT));
	cylinder->setColour(RGB(0, 0, 1));
	cylinder->setDiffuse(RGB(0.2));
	cylinder->shift(650, 2000, 400);
	cylinder->setPhong(20.0);
	
	std::cout << "Cylinder created\n";
	
	cone->stretch(1, 1, 2);
	coneT->rotate(X_AXIS, 180);
	coneB->shift(0, 0, -600);
	CSGPtr coneOb;
	coneOb = combine(std::move(cone), SUB, std::move(coneT));
	coneOb = combine(std::move(coneOb), SUB, std::move(coneB));
	coneOb->setColour(RGB(0, 1, 0));
	coneOb->setDiffuse(RGB(0.4));
	coneOb->setOpacity(RGB(0.9), 1.5);
	coneOb->setPhong(50);
	coneOb->shift(-650, 2000, 700);
	
	std::cout << "Cone done\n";
	
	CSGPtr tree;
	tree = combine(std::move(basePlane), ADD, std::move(sp));
	tree = combine(std::move(tree), ADD, std::move(cylinder));
	//tree = combine(std::move(tree), ADD, std::move(coneOb));
	tree = combine(std::move(tree), ADD, std::move(cube));
	
	printCSGTree(tree.get());

	
	s.objectTree = std::move(tree);
	
	s.ambient = RGB(0.2);
	RGB lightInt(0.6);
	Vector lightDir(-1, 1, -1);
	lightDir = lightDir.normalised();
	s.addMoon(lightInt, lightDir);
	lightDir = Vector(1, 1, -0.5);
	lightDir = lightDir.normalised();
	s.addMoon(lightInt, lightDir);
	
	s.viewPoint = Point(1000, 500, 500);
	s.toPoint = Point(0, 2000, 500);
	
	//s.screen = Screen(320, 180);
	s.screen = Screen(1280, 720);
	
	s.transmissionCount = 5;
	
	//render(s, "testCube");
	
	for (int frame = 0; frame < nFrames; frame++) {
	
		frame = 25;

		double xpos, ypos, zpos;
		//xpos = 0.0 - 4000.0 * std::sin(2 * M_PI * frame / 100.0);
		xpos = 0.0 - 3000.0 * std::sin(2 * M_PI * frame / 100.0);
		ypos = 2000.0 - 2000.0 * std::cos(2 * M_PI * frame / 100.0);
		//zpos = 500.0 + 1000.0 * (1 - std::cos(2 * M_PI * frame / 100.0));
		zpos = 500.0;
		s.viewPoint.coord[0] = xpos;
		s.viewPoint.coord[1] = ypos;
		s.viewPoint.coord[2] = zpos;
		std::string fname = std::to_string(frame) + "L";
		render (s, std::string(4 - fname.length(), '0') + fname);
		
		//frame = 101;
	}
	
	return 0;
}



