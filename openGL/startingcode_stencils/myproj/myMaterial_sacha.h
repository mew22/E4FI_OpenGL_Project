#pragma once
#include <GL/glew.h>
#include <string>

class myMaterial
{
public:
	GLfloat material_Ka[4];
	GLfloat material_Kd[4];
	GLfloat material_Ks[4];
	GLfloat material_Sh;
	std::string texture_name;
	std::string material_name;

	void setMaterial(GLfloat ka_r, GLfloat ka_g, GLfloat ka_b, GLfloat ka_a,
		GLfloat kd_r, GLfloat kd_g, GLfloat kd_b, GLfloat kd_a,
		GLfloat ks_r, GLfloat ks_g, GLfloat ks_b, GLfloat ks_a,
		GLfloat s);

	void setMaterial(const myMaterial & m);

	void activateMaterial();

	myMaterial(void);
	~myMaterial(void);
};

