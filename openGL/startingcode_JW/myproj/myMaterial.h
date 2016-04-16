class myMaterial
{
	public:
	GLfloat material_Ka[4]; //ambient color
	GLfloat material_Kd[4]; //diffuse color
	GLfloat material_Ks[4]; //specular color
	GLfloat material_Sh; //shininess coefficient
	string material_name; //name of the material
};