#pragma once
#include <GL/glew.h>
#include <string>

class myTexture
{
public:
	int width, height, pixelsize;
	bool repeat = false;
	float scaleU = 1, scaleV = 1;
	GLuint texName;
	std::string strName;


	myTexture();
	GLubyte * readFile(char *filename, int &, int &);
	bool readTexture(char *filename);
	bool readCubemap(char *, char *, char *, char *, char *, char *);

	void activateTexture(GLenum texture_slot, GLenum texture_type = GL_TEXTURE_2D);
};

