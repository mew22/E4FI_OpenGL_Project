#pragma once

#include <math.h>
#include <GL/glew.h>
#include <vector>
#include <string>
#include <fstream>
#include "vector3d.h"
#include "myMaterial.h"
#include "myTexture.h"

class mySubObject3D
{
public:
	myTexture *texture = nullptr, *bump = nullptr, *cubemap = nullptr;
	myMaterial *material = nullptr;
	int start, end; //then indices are (start, end-1).
	string name;

	mySubObject3D(int s, int e) { start = s; end = e; }
};