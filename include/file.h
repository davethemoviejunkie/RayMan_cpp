/***************************************************************************

file.h

Header file for writing files

Rewrite of Python Raytracer in C++

***************************************************************************/

#ifndef FILE_H
#define FILE_H

#include <vector>

#include "colour.h"

void outputPicture(std::vector<RGB>, int, int, std::string);

#endif
