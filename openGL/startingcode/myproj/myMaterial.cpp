#include "StdAfx.h"
#include "myMaterial.h"


myMaterial::myMaterial(void)
{
	setMaterial(0.1,0.1,0.1,1,	1,0,1,1,	0.2,0.2,0.2,1,	40);
}


myMaterial::~myMaterial(void)
{
}

void myMaterial::setMaterial(GLfloat ka_r, GLfloat ka_g, GLfloat ka_b, GLfloat ka_a, 
	GLfloat kd_r, GLfloat kd_g, GLfloat kd_b, GLfloat kd_a, 
	GLfloat ks_r, GLfloat ks_g, GLfloat ks_b, GLfloat ks_a,
	GLfloat s)
{
	material_Ka[0] = ka_r; material_Ka[1] = ka_g; material_Ka[2] = ka_b; material_Ka[3] = ka_a;
	material_Kd[0] = kd_r; material_Kd[1] = kd_g; material_Kd[2] = kd_b; material_Kd[3] = kd_a;
	material_Ks[0] = ks_r; material_Ks[1] = ks_g; material_Ks[2] = ks_b; material_Ks[3] = ks_a;
	material_Sh = s;
}

void myMaterial::setMaterial(const myMaterial & m)
{
	setMaterial(m.material_Ka[0], m.material_Ka[1], m.material_Ka[2], m.material_Ka[3],
				m.material_Kd[0], m.material_Kd[1], m.material_Kd[2], m.material_Kd[3], 
				m.material_Ks[0], m.material_Ks[1], m.material_Ks[2], m.material_Ks[3], m.material_Sh);
}

void myMaterial::activateMaterial()
{
	glMaterialfv(GL_FRONT, GL_AMBIENT, material_Ka);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, material_Kd);
	glMaterialfv(GL_FRONT, GL_SPECULAR, material_Ks);
	glMaterialf(GL_FRONT, GL_SHININESS, material_Sh);
}
