#pragma once
#include <GL/glew.h>

class myTexture
{
public:
	int width, height, pixelsize;
	GLuint texName;

	myTexture();
	GLubyte * readFile(char *filename, int &, int &);
	bool readTexture(char *filename);
	bool readCubemap(char *, char *, char *, char *, char *, char *);

	void activateTexture(GLenum texture_slot, GLenum texture_type = GL_TEXTURE_2D);
};

