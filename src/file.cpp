/***************************************************************************

file.cpp

File routines -- outputting the picture

Rewrite of Python Raytracer in C++

***************************************************************************/

#include <fstream>

#include "file.h"

void outputPicture(std::vector<RGB> theImg, int width, int height, std::string filename) {

	std::ofstream imageFile;
	
	int across, down, thisPix;

	imageFile.open(filename + ".ppm");

	imageFile << "P3\n";
	imageFile << width << " " << height << "\n"; // dimensions
	imageFile << "255\n"; // Max pixel

	for (int down = 0; down < height; down++)
		for (int across = 0; across < width; across++) {
			thisPix = down * width + across;
			for (int i = 0; i < COLOURSPACE; i++)
				imageFile << (int)(theImg[thisPix].value[i] * 255 + 0.5) << " ";
			imageFile << "\n";
		}
	imageFile.close();

};
