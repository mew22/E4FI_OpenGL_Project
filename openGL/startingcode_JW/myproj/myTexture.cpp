#include "StdAfx.h"
#include <fstream>
#include <iostream>
#include "myTexture.h"

myTexture::myTexture()
{
	texName = -1;
}

GLubyte * myTexture::readFile(char *filename, int & w, int & h)
{
	FILE *inFile; 
	char buffer[100]; 
    GLubyte *mytexture; 
	unsigned char c; 
	int maxVal;

	if( (inFile = fopen(filename, "rb")) == NULL) {
		std::cout << "Unable to open the ppm file";
		getchar();
		exit(1);
	}


	//Read file type identifier (magic number)
	fgets(buffer, sizeof(buffer), inFile);
	if ((buffer[0] != 'P') || (buffer[1] != '6')) {
		fprintf (stderr, "not a binary ppm file %s\n", filename);
		return 0;
	}

	if(buffer[2] == 'A')
		pixelsize = 4;
	else
		pixelsize = 3;

	//Read image size
	do fgets(buffer, sizeof (buffer), inFile);
	while (buffer[0] == '#');
	sscanf (buffer, "%d %d", &w, &h);

	//Read maximum pixel value (usually 255)
	do fgets (buffer, sizeof (buffer), inFile);
	while (buffer[0] == '#');
	sscanf (buffer, "%d", &maxVal);

	//Allocate RGBA texture buffer
	int memSize = w * h * 4 * sizeof(GLubyte);
	mytexture = new GLubyte[memSize];

	// read RGB data and set alpha value
	for (int i = 0; i < memSize; i++) {
		if ((i % 4) < 3 || pixelsize == 4) {
			c = fgetc(inFile);
			mytexture[i]=(GLubyte) c;
        }
		else mytexture[i] = (GLubyte) 255; //Set alpha to opaque
    }
    fclose(inFile);
	return mytexture;
}

bool myTexture::readTexture(char *filename)
{ 
	GLubyte *mytexture = readFile(filename, width, height);
	glGenTextures(1, &texName) ; 
    glBindTexture (GL_TEXTURE_2D, texName) ; 

    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR) ; 
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR) ; 
   
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, (GLuint)width, (GLuint)height, 0, GL_RGBA, GL_UNSIGNED_BYTE, mytexture);

	delete[] mytexture;
	return true;
}


bool myTexture::readCubemap(char *filenamepx, char *filenamenx, char *filenamepy, char *filenameny, char *filenamepz, char *filenamenz)
{ 
	GLubyte *cubemap;
	int w, h;

	glGenTextures(1, &texName) ; 
    glBindTexture (GL_TEXTURE_CUBE_MAP, texName) ; 

	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	
	cubemap = readFile(filenamepx, w, h);
	glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X, 0, GL_RGBA, (GLuint)w, (GLuint)h,0, GL_RGBA, GL_UNSIGNED_BYTE, cubemap);
	delete[] cubemap;

	cubemap = readFile(filenamenx, w, h);
	glTexImage2D(GL_TEXTURE_CUBE_MAP_NEGATIVE_X, 0, GL_RGBA, (GLuint)w, (GLuint)h,0, GL_RGBA, GL_UNSIGNED_BYTE, cubemap);
	delete[] cubemap;

	cubemap = readFile(filenamepy, w, h);
	glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_Y, 0, GL_RGBA, (GLuint)w, (GLuint)h,0, GL_RGBA, GL_UNSIGNED_BYTE, cubemap);
	delete[] cubemap;

	cubemap = readFile(filenameny, w, h);
	glTexImage2D(GL_TEXTURE_CUBE_MAP_NEGATIVE_Y, 0, GL_RGBA, (GLuint)w, (GLuint)h,0, GL_RGBA, GL_UNSIGNED_BYTE, cubemap);
	delete[] cubemap;

	cubemap = readFile(filenamepz, w, h);
	glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_Z, 0, GL_RGBA, (GLuint)w, (GLuint)h,0, GL_RGBA, GL_UNSIGNED_BYTE, cubemap);
	delete[] cubemap;

	cubemap = readFile(filenamenz, w, h);
	glTexImage2D(GL_TEXTURE_CUBE_MAP_NEGATIVE_Z, 0, GL_RGBA, (GLuint)w, (GLuint)h,0, GL_RGBA, GL_UNSIGNED_BYTE, cubemap);
	delete[] cubemap;

	width = w; height = h;
	return true;
}


void myTexture::activateTexture(GLenum texture_slot, GLenum texture_type)
{
	glActiveTexture(texture_slot); 
	glBindTexture(texture_type, texName) ; 
}

