#pragma once
#include <math.h>
#include <GL/glew.h>
#include "myTexture.h"
#include "myMaterial.h"
#include <vector>
#include <string>
#include "vector3d.h"
class mySubObject3D
{
public:
	myMaterial material;
	int start_index, end_index;
	myTexture texture, bump, cubemap;
	int start, end; //then indices are (start, end-1).
	string name;

	mySubObject3D(int s, int e) { start = s; end = e; }
};