#pragma once
#include <math.h>
#include <GL/glew.h>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include "vector3d.h"
#include "myMaterial.h"
#include "mySubObject3D.h"

#define PI 3.14159265

using namespace std;



class myObject3D
{
public:
	GLuint buffers[5];

	GLuint myshaderprogram; //the ID of shader program, in case multiple shaders
	myMaterial *material; //contains the color description of the object.
	vector<GLfloat> vertices;
	vector<GLuint> indices;
	vector<GLfloat> normals;
	vector<GLfloat> texturecoordinates;
	vector<myTexture*> texTmp;
	vector<GLfloat> tangents;
	map<string, myMaterial *> materials;
	

	vector<mySubObject3D *> parts; //contains subparts of the scene.

	//Model matrix
	glm::mat4 model_matrix;
	
	myObject3D(GLuint shader) {
		clear();
		myshaderprogram = shader;
		model_matrix = glm::mat4(1.0f);

		material = (myMaterial*)malloc(sizeof(myMaterial));
		materials = map<string, myMaterial*>();

	}

	void clear() {
		vertices.clear();
		indices.clear();
		normals.clear();
		texturecoordinates.clear();
		parts.clear();
		tangents.clear();
	}
 
	void readMesh(char *filename)
	{
		string s, t;
		string tmp;
		ifstream fin(filename);

		if(!fin.is_open()) cout << "Unable to open the file";

		while (getline(fin, s))
		{
			stringstream myline(s);
			myline >> t;
			if(t == "v")
			{
				myline >> tmp;
				double x = stof(tmp.substr(0, tmp.find("/")));
				//std::cout << "(" << stof(tmp.substr(0, tmp.find("/")));

				myline >> tmp;
				double y = stof(tmp.substr(0, tmp.find("/")));
				//std::cout << ", " <<  stof(tmp.substr(0, tmp.find("/")));

				myline >> tmp;
				double z = stof(tmp.substr(0, tmp.find("/")));
				//std::cout << ", " <<  stof(tmp.substr(0, tmp.find("/"))) << ")\n";

				vertices.push_back(x);
				vertices.push_back(y);
				vertices.push_back(z);
			}

			if(t == "f")
			{
				myline >> tmp;
				int f1 = stoi(tmp.substr(0, tmp.find("/")));

				myline >> tmp;
				int f2 = stoi(tmp.substr(0, tmp.find("/")));

				myline >> tmp;
				int f3 = stoi(tmp.substr(0, tmp.find("/")));

				indices.push_back(f1-1);
				indices.push_back(f2-1);
				indices.push_back(f3-1);
			}
		}
		normalize();
	}

	void normalize()
	{
		int i;
		int tmpxmin = 0, tmpymin = 0, tmpzmin = 0, tmpxmax = 0, tmpymax = 0, tmpzmax = 0;

		int n = vertices.size()/3;

		for (i=0;i<n;i++) {
			if (vertices[3*i] < vertices[3*tmpxmin]) tmpxmin = i;
			if (vertices[3*i] > vertices[3*tmpxmax]) tmpxmax = i;

			if (vertices[3*i+1] < vertices[3*tmpymin+1]) tmpymin = i;
			if (vertices[3*i+1] > vertices[3*tmpymax+1]) tmpymax = i;

			if (vertices[3*i+2] < vertices[3*tmpzmin+2]) tmpzmin = i;
			if (vertices[3*i+2] > vertices[3*tmpzmax+2]) tmpzmax = i;
		}

		double xmin = vertices[3*tmpxmin], xmax = vertices[3*tmpxmax], 
			   ymin = vertices[3*tmpymin+1], ymax = vertices[3*tmpymax+1], 
			   zmin = vertices[3*tmpzmin+2], zmax = vertices[3*tmpzmax+2];

		double scale = (xmax-xmin) <= (ymax-ymin) ? (xmax-xmin) : (ymax-ymin);
		//double scale = fmin( (xmax-xmin), (ymax-ymin) );
		scale = scale >= (zmax-zmin) ? scale : (zmax-zmin);
		//scale = fmax(scale, (zmax-zmin));

		for (i=0;i<n;i++) {
			vertices[3*i] -= (xmax+xmin)/2;
			vertices[3*i+1] -= (ymax+ymin)/2;
			vertices[3*i+2] -= (zmax+zmin)/2;

			vertices[3*i] /= scale;
			vertices[3*i+1] /= scale;
			vertices[3*i+2] /= scale;
		}
	}

	void computeNormal(int v1, int v2, int v3, float & x, float & y, float & z)
	{
		double dx1 = vertices[v2*3] - vertices[v1*3];
		double dx2 = vertices[v3*3] - vertices[v2*3];
		double dy1 = vertices[v2*3+1] - vertices[v1*3+1];
		double dy2 = vertices[v3*3+1] - vertices[v2*3+1];
		double dz1 = vertices[v2*3+2] - vertices[v1*3+2];
		double dz2 = vertices[v3*3+2] - vertices[v2*3+2];


		double dx = dy1 * dz2 - dz1 * dy2;
		double dy = dz1 * dx2 - dx1 * dz2;
		double dz = dx1 * dy2 - dy1 * dx2;

		double length = sqrt(dx*dx + dy*dy + dz*dz);
		if (length <= 0)
		{
			//cout << "Error! vector length is zero\n";
			x = y = z = 1.0f;
			return;
		}

		x = dx/length;
		y = dy/length;
		z = dz/length;
	}

	void computeNormals( )
	{
		normals.clear();
		int i, j;
		float x1, y1, z1;

		int n = vertices.size()/3;
		int m = indices.size()/3;

		normals.resize(3*n);
		int *incidences = new int[n];
		for (i=0;i<3*n;i++) normals[i] = 0.0;
		for (i=0;i<n;i++) incidences[i] = 0;

		for (j=0;j<m;j++)
		{
			computeNormal(indices[3*j], indices[3*j+1], indices[3*j+2], x1, y1, z1);
			normals[3*indices[3*j]] += x1; normals[3*indices[3*j]+1] += y1; normals[3*indices[3*j]+2] += z1;
			normals[3*indices[3*j+1]] += x1; normals[3*indices[3*j+1]+1] += y1; normals[3*indices[3*j+1]+2] += z1;
			normals[3*indices[3*j+2]] += x1; normals[3*indices[3*j+2]+1] += y1; normals[3*indices[3*j+2]+2] += z1;
			incidences[indices[3*j]]++; incidences[indices[3*j+1]]++; incidences[indices[3*j+2]]++;
		}
		for (i=0;i<n;i++)
			if (incidences[i]!=0) 
			{
				normals[3*i] /= incidences[i]; normals[3*i+1] /= incidences[i]; normals[3*i+2] /= incidences[i];
			}
	}


	void displayObject(GLuint shaderprogram, glm::mat4 viewmatrix)
	{
		glUniformMatrix4fv(glGetUniformLocation(shaderprogram, "mymodel_matrix"), 1, GL_FALSE, &model_matrix[0][0]);
		glm::mat3 normal_matrix = glm::transpose(glm::inverse(glm::mat3(viewmatrix*model_matrix)));
		glUniformMatrix3fv(glGetUniformLocation(shaderprogram, "mynormal_matrix"), 1, GL_FALSE, &normal_matrix[0][0]);

		glUniform4fv(glGetUniformLocation(shaderprogram, "mymodel_kd"), 1, material->material_Kd);
		glUniform4fv(glGetUniformLocation(shaderprogram, "mymodel_ka"), 1, material->material_Ka);
		glUniform4fv(glGetUniformLocation(shaderprogram, "mymodel_ks"), 1,material->material_Ks);
		glUniform1f(glGetUniformLocation(shaderprogram, "mymodel_ns"), material->material_Sh);


		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

		glEnableVertexAttribArray(2);
		glBindBuffer(GL_ARRAY_BUFFER, buffers[2]);
		glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, 0);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffers[1]);
		glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);

		glUniformMatrix4fv(glGetUniformLocation(shaderprogram, "mymodel_matrix"),
		1, GL_FALSE, &model_matrix[0][0]);
	}

	void displayNormals()
	{
		glBegin(GL_LINES);
		for (int i = 0; i < this->normals.size(); i += 3)
		{
			glVertex3f(vertices[i], vertices[i + 1], vertices[i + 2]);
			glVertex3f(vertices[i] + 0.1*normals[i],
			vertices[i + 1] + 0.1*normals[i + 1],
			vertices[i + 2] + 0.1*normals[i + 2]);
		}
		glEnd();
	}

	void translate(double x, double y, double z)
	{
		glm::mat4 tmp = glm::translate(glm::vec3(x,y,z));
		model_matrix = tmp * model_matrix;
	}
	void rotate(double axis_x, double axis_y, double axis_z, double angle)
	{
		glm::mat4 tmp = glm::rotate((float) angle, glm::vec3(axis_x, axis_y, axis_z));
		model_matrix = tmp * model_matrix;
	}
	void scale(double axis_x, double axis_y, double axis_z)
	{
		glm::mat4 tmp = glm::scale(glm::vec3(axis_x, axis_y, axis_z));
		model_matrix = tmp * model_matrix;
	}

	void normalize(vector<GLfloat> & myvertices)
	{
		if (myvertices.size() < 4) return;
		int i;
		int tmpxmin = 0, tmpymin = 0, tmpzmin = 0, tmpxmax = 0, tmpymax = 0, tmpzmax = 0;

		int n = myvertices.size() / 3;

		for (i = 0; i<n; i++) {
			if (myvertices[3 * i] < myvertices[3 * tmpxmin]) tmpxmin = i;
			if (myvertices[3 * i] > myvertices[3 * tmpxmax]) tmpxmax = i;

			if (myvertices[3 * i + 1] < myvertices[3 * tmpymin + 1]) tmpymin = i;
			if (myvertices[3 * i + 1] > myvertices[3 * tmpymax + 1]) tmpymax = i;

			if (myvertices[3 * i + 2] < myvertices[3 * tmpzmin + 2]) tmpzmin = i;
			if (myvertices[3 * i + 2] > myvertices[3 * tmpzmax + 2]) tmpzmax = i;
		}

		double xmin = myvertices[3 * tmpxmin], xmax = myvertices[3 * tmpxmax],
			ymin = myvertices[3 * tmpymin + 1], ymax = myvertices[3 * tmpymax + 1],
			zmin = myvertices[3 * tmpzmin + 2], zmax = myvertices[3 * tmpzmax + 2];

		double scale = (xmax - xmin) <= (ymax - ymin) ? (xmax - xmin) : (ymax - ymin);
		scale = scale >= (zmax - zmin) ? scale : (zmax - zmin);

		for (i = 0; i<n; i++) {
			myvertices[3 * i] -= (xmax + xmin) / 2;
			myvertices[3 * i + 1] -= (ymax + ymin) / 2;
			myvertices[3 * i + 2] -= (zmax + zmin) / 2;

			myvertices[3 * i] /= scale;
			myvertices[3 * i + 1] /= scale;
			myvertices[3 * i + 2] /= scale;
		}
	}
	void computeTangent(int v0, int v1, int v2, float & x, float & y, float & z)
	{
		float du1 = texturecoordinates[2 * v1] - texturecoordinates[2 * v0];
		float dv1 = texturecoordinates[2 * v1 + 1] - texturecoordinates[2 * v0 + 1];
		float du2 = texturecoordinates[2 * v2] - texturecoordinates[2 * v0];
		float dv2 = texturecoordinates[2 * v2 + 1] - texturecoordinates[2 * v0 + 1];

		float f = 1.0f / (du1 * dv2 - du2 * dv1);
		if ((du1*dv2 - du2*dv1) == 0){
			x = y = z = 0; return;
		}

		float e1x = vertices[3 * v1] - vertices[3 * v0];
		float e1y = vertices[3 * v1 + 1] - vertices[3 * v0 + 1];
		float e1z = vertices[3 * v1 + 2] - vertices[3 * v0 + 2];

		float e2x = vertices[3 * v2] - vertices[3 * v0];
		float e2y = vertices[3 * v2 + 1] - vertices[3 * v0 + 1];
		float e2z = vertices[3 * v2 + 2] - vertices[3 * v0 + 2];

		x = f * (dv2 * e1x - dv1 * e2x);
		y = f * (dv2 * e1y - dv1 * e2y);
		z = f * (dv2 * e1z - dv1 * e2z);
	}

	void computeTangents()
	{
		int i, j, k;
		GLfloat x1, y1, z1;

		int n = vertices.size() / 3;
		int m = indices.size() / 3;

		tangents.resize(3 * n);
		int *incidences = new int[n];
		for (i = 0; i<3 * n; i++) tangents[i] = 0.0;
		for (i = 0; i<n; i++) incidences[i] = 0;

		for (j = 0; j<m; j++)
		{
			computeTangent(indices[3 * j], indices[3 * j + 1], indices[3 * j + 2], x1, y1, z1);
			tangents[3 * indices[3 * j]] += x1; tangents[3 * indices[3 * j] + 1] += y1; tangents[3 * indices[3 * j] + 2] += z1;
			tangents[3 * indices[3 * j + 1]] += x1; tangents[3 * indices[3 * j + 1] + 1] += y1; tangents[3 * indices[3 * j + 1] + 2] += z1;
			tangents[3 * indices[3 * j + 2]] += x1; tangents[3 * indices[3 * j + 2] + 1] += y1; tangents[3 * indices[3 * j + 2] + 2] += z1;
			incidences[indices[3 * j]]++; incidences[indices[3 * j + 1]]++; incidences[indices[3 * j + 2]]++;
		}
		for (i = 0; i<n; i++) {
			float l = sqrt(tangents[3 * i] * tangents[3 * i] + tangents[3 * i + 1] * tangents[3 * i + 1] + tangents[3 * i + 2] * tangents[3 * i + 2]);
			tangents[3 * i] /= l; tangents[3 * i + 1] /= l; tangents[3 * i + 2] /= l;
		}
	}
	void readScene(char *filename) {
		string s, t;
		string tmp;
		ifstream fin(filename);

		if (!fin.is_open()) cout << "Unable to open the file";

		mySubObject3D* currentPart = nullptr;
		std::vector<GLfloat> tmp_normals;
		std::vector<GLfloat> tmp_textures;

		string part_name = "";

		while (getline(fin, s))
		{
			stringstream myline(s);
			myline >> t;
			if (t == "g") {
				part_name = "";
				while (myline >> tmp) {
					if (part_name != "")
						part_name.append(" ");
					part_name.append(tmp);
				}

			}

			if (t.compare("usemtl") == 0) {

				// create part
				currentPart = new mySubObject3D();
				int end_index = indices.size() - 1;
				currentPart->start = indices.size();
				if (parts.size() > 0) {
					parts.at(parts.size() - 1)->end = end_index;
					// remove old if no faces
					if (parts.at(parts.size() - 1)->end <= parts.at(parts.size() - 1)->start)
						parts.pop_back();
				}
				parts.push_back(currentPart);
				currentPart->object_name = part_name;

				myline >> tmp;
				if (materials[tmp] != nullptr) {
					currentPart->tex = new myTexture();
					if(materials[tmp]->texture_name != "")
						currentPart->tex->readTexture((char*)materials[tmp]->texture_name.c_str());
					currentPart->material = materials[tmp];
				}

				//currentPart->bump_tex->readTexture("Ppm/4351-normal.ppm");
			}

			if (t.compare("mtllib") == 0) {
				myline >> tmp;
				readMaterials(tmp);
			}


			if (t == "v")
			{
				myline >> tmp;
				double x = stof(tmp.substr(0, tmp.find("/")));
				//std::cout << "(" << stof(tmp.substr(0, tmp.find("/")));

				myline >> tmp;
				double y = stof(tmp.substr(0, tmp.find("/")));
				//std::cout << ", " <<  stof(tmp.substr(0, tmp.find("/")));

				myline >> tmp;
				double z = stof(tmp.substr(0, tmp.find("/")));
				//std::cout << ", " <<  stof(tmp.substr(0, tmp.find("/"))) << ")\n";

				vertices.push_back(x);
				vertices.push_back(y);
				vertices.push_back(z);
			}

			if (t == "vn")
			{
				myline >> tmp;
				double x = stof(tmp.substr(0, tmp.find("/")));
				//std::cout << "(" << stof(tmp.substr(0, tmp.find("/")));

				myline >> tmp;
				double y = stof(tmp.substr(0, tmp.find("/")));
				//std::cout << ", " <<  stof(tmp.substr(0, tmp.find("/")));

				myline >> tmp;
				double z = stof(tmp.substr(0, tmp.find("/")));
				//std::cout << ", " <<  stof(tmp.substr(0, tmp.find("/"))) << ")\n";

				normals.push_back(x);
				normals.push_back(y);
				normals.push_back(z);
			}

			if (t == "vt")
			{
				myline >> tmp;
				double x = stof(tmp.substr(0, tmp.find("/")));
				//std::cout << "(" << stof(tmp.substr(0, tmp.find("/")));

				myline >> tmp;
				double y = stof(tmp.substr(0, tmp.find("/")));
				//std::cout << ", " <<  stof(tmp.substr(0, tmp.find("/")));

				texturecoordinates.push_back(x);
				texturecoordinates.push_back(y);
			}

			if (t == "f")
			{

				int nb_vertices = vertices.size() / 3;
				tmp_normals.resize(nb_vertices * 3);
				tmp_textures.resize(nb_vertices * 2);

				myline >> tmp;
				int vx[3] = { 0, 0, 0 };
				vx[0] = atoi((tmp.substr(0, tmp.find("/"))).c_str());
				if (tmp.find("/") != -1)
					vx[1] = atoi((tmp.substr(tmp.find("/") + 1, tmp.rfind("/"))).c_str());
				if (tmp.rfind("/") != -1 && tmp.rfind("/") != tmp.find("/"))
					vx[2] = atoi((tmp.substr(tmp.rfind("/") + 1, tmp.size())).c_str());

				int f = vx[0];

				if (vx[1] != 0) {
					tmp_textures.at((vx[0] - 1) * 2) = texturecoordinates.at((vx[1] - 1) * 2);
					tmp_textures.at((vx[0] - 1) * 2 + 1) = texturecoordinates.at((vx[1] - 1) * 2 + 1);
				}

				if (vx[2] != 0) {
					tmp_normals.at((vx[0] - 1) * 3) = normals.at((vx[2] - 1) * 3);
					tmp_normals.at((vx[0] - 1) * 3 + 1) = normals.at((vx[2] - 1) * 3 + 1);
					tmp_normals.at((vx[0] - 1) * 3 + 2) = normals.at((vx[2] - 1) * 3 + 2);
				}

				myline >> tmp;
				vx[0] = atoi((tmp.substr(0, tmp.find("/"))).c_str());
				if (tmp.find("/") != -1)
					vx[1] = atoi((tmp.substr(tmp.find("/") + 1, tmp.rfind("/"))).c_str());
				if (tmp.rfind("/") != -1 && tmp.rfind("/") != tmp.find("/"))
					vx[2] = atoi((tmp.substr(tmp.rfind("/") + 1, tmp.size())).c_str());
				int l = vx[0];

				if (vx[1] != 0) {
					tmp_textures.at((vx[0] - 1) * 2) = texturecoordinates.at((vx[1] - 1) * 2);
					tmp_textures.at((vx[0] - 1) * 2 + 1) = texturecoordinates.at((vx[1] - 1) * 2 + 1);
				}

				if (vx[2] != 0) {
					tmp_normals.at((vx[0] - 1) * 3) = normals.at((vx[2] - 1) * 3);
					tmp_normals.at((vx[0] - 1) * 3 + 1) = normals.at((vx[2] - 1) * 3 + 1);
					tmp_normals.at((vx[0] - 1) * 3 + 2) = normals.at((vx[2] - 1) * 3 + 2);
				}


				while (myline >> tmp)
				{
					indices.push_back(f - 1);
					indices.push_back(l - 1);

					vx[0] = atoi((tmp.substr(0, tmp.find("/"))).c_str());
					if (tmp.find("/") != -1)
						vx[1] = atoi((tmp.substr(tmp.find("/") + 1, tmp.rfind("/"))).c_str());
					if (tmp.rfind("/") != -1 && tmp.rfind("/") != tmp.find("/"))
						vx[2] = atoi((tmp.substr(tmp.rfind("/") + 1, tmp.size())).c_str());
					l = vx[0];

					if (vx[1] != 0) {
						tmp_textures.at((vx[0] - 1) * 2) = texturecoordinates.at((vx[1] - 1) * 2);
						tmp_textures.at((vx[0] - 1) * 2 + 1) = texturecoordinates.at((vx[1] - 1) * 2 + 1);
					}

					if (vx[2] != 0) {
						tmp_normals.at((vx[0] - 1) * 3) = normals.at((vx[2] - 1) * 3);
						tmp_normals.at((vx[0] - 1) * 3 + 1) = normals.at((vx[2] - 1) * 3 + 1);
						tmp_normals.at((vx[0] - 1) * 3 + 2) = normals.at((vx[2] - 1) * 3 + 2);
					}

					indices.push_back(l - 1);
				}
			}
		}
		parts.at(parts.size() - 1)->end= indices.size() - 1;
		for (int i = 0; i < parts.size(); ++i) {
			if (parts.at(i)->material == nullptr)
				parts.at(i)->material = new myMaterial();
		}
		computeNormals();
		//normals = std::vector<GLfloat>(tmp_normals);
		texturecoordinates = std::vector<GLfloat>(tmp_textures);
	}

	void readMaterials(string filename) {
		string s, t;
		string tmp;
		string matName;
		ifstream fin(filename);

		if (!fin.is_open()) cout << "Unable to open the file";

		while (getline(fin, s))
		{
			stringstream myline(s);
			myline >> t;

			if (t.compare("newmtl") == 0) {
				myline >> matName;
				materials[matName] = new myMaterial();
				materials[matName]->material_name = matName;
				continue;
			}

			if (t.compare("Kd") == 0) {
				myline >> tmp;
				materials[matName]->material_Kd[0] = stof(tmp);
				myline >> tmp;
				materials[matName]->material_Kd[1] = stof(tmp);
				myline >> tmp;
				materials[matName]->material_Kd[2] = stof(tmp);
				materials[matName]->material_Kd[3] = 1;
			}

			if (t.compare("map_Kd") == 0) {
				myline >> tmp;
				materials[matName]->texture_name.append(tmp);
			}

			if (t.compare("Ks") == 0) {
				myline >> tmp;
				materials[matName]->material_Ks[0] = stof(tmp);
				myline >> tmp;
				materials[matName]->material_Ks[1] = stof(tmp);
				myline >> tmp;
				materials[matName]->material_Ks[2] = stof(tmp);
				materials[matName]->material_Ks[3] = 1;
			}

			if (t.compare("Ka") == 0) {
				myline >> tmp;
				materials[matName]->material_Ka[0] = stof(tmp);
				myline >> tmp;
				materials[matName]->material_Ka[1] = stof(tmp);
				myline >> tmp;
				materials[matName]->material_Ka[2] = stof(tmp);
				materials[matName]->material_Ka[3] = 1;
			}

			if (t.compare("Ns") == 0) {
				myline >> tmp;
				materials[matName]->material_Sh = stof(tmp);
			}
		}
	}

	void createObjectBuffers()
	{
		glGenBuffers(4, buffers);

		glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
		glBufferData(GL_ARRAY_BUFFER, vertices.size() * 4, &vertices.front(), GL_STATIC_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, buffers[1]);
		glBufferData(GL_ARRAY_BUFFER, indices.size() * 4, &indices.front(), GL_STATIC_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, buffers[2]);
		glBufferData(GL_ARRAY_BUFFER, normals.size() * 4, &normals.front(), GL_STATIC_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, buffers[3]);
		glBufferData(GL_ARRAY_BUFFER, texturecoordinates.size() * 4, &texturecoordinates.front(), GL_STATIC_DRAW);


		//glBindBuffer(GL_ARRAY_BUFFER, buffers[4]);
		//glBufferData(GL_ARRAY_BUFFER, tangents.size() * 4, &tangents.front(), GL_STATIC_DRAW);
	}

	void displayScene(glm::mat4 viewmatrix){
		glUniformMatrix4fv(glGetUniformLocation(myshaderprogram, "mymodel_matrix"), 1, GL_FALSE, &model_matrix[0][0]);
		glm::mat3 normal_matrix = glm::transpose(glm::inverse(glm::mat3(viewmatrix*model_matrix)));
		glUniformMatrix3fv(glGetUniformLocation(myshaderprogram, "mynormal_matrix"), 1, GL_FALSE, &normal_matrix[0][0]);

		glEnableVertexAttribArray(0);
		glEnableVertexAttribArray(2);
		glEnableVertexAttribArray(3);
		//glEnableVertexAttribArray(4);

		glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

		
		glBindBuffer(GL_ARRAY_BUFFER, buffers[2]);
		glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, 0);

		glBindBuffer(GL_ARRAY_BUFFER, buffers[3]);
		glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, 0, 0);

		//glBindBuffer(GL_ARRAY_BUFFER, buffers[4]);
		//glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, 0, 0);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffers[1]);
		//glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);

		for (int i = 0; i < parts.size(); ++i) {

			mySubObject3D sub_obj = *(parts.at(i));

			glUniform4fv(glGetUniformLocation(myshaderprogram, "Ks"), 1, &(sub_obj.material->material_Ks[0]));
			glUniform4fv(glGetUniformLocation(myshaderprogram, "Kd"), 1, &(sub_obj.material->material_Kd[0]));
			glUniform4fv(glGetUniformLocation(myshaderprogram, "Ka"), 1, &(sub_obj.material->material_Ka[0]));
			glUniform1f(glGetUniformLocation(myshaderprogram, "Ns"), sub_obj.material->material_Sh);
			glUniform1i(glGetUniformLocation(myshaderprogram, "wave"), 0);
			glUniform1i(glGetUniformLocation(myshaderprogram, "using_textures"), 0);
			if (sub_obj.object_name == "wave scene")
				glUniform1i(glGetUniformLocation(myshaderprogram, "wave"), 1);

			if (sub_obj.material->texture_name != "") {
				glUniform1i(glGetUniformLocation(myshaderprogram, "using_textures"), 1);

				glActiveTexture(GL_TEXTURE8);
				glBindTexture(GL_TEXTURE_2D, sub_obj.tex->texName);
				glUniform1i(glGetUniformLocation(myshaderprogram, "tex"), 8);

				//glActiveTexture(GL_TEXTURE9);
				//glBindTexture(GL_TEXTURE_2D, sub_obj.bump_tex->texName);
				//glUniform1i(glGetUniformLocation(shaderprogram,"bump_tex"),9);
			}

			glDrawElements(GL_TRIANGLES, sub_obj.end - sub_obj.start + 1, GL_UNSIGNED_INT, (void*)(sub_obj.start * sizeof(GLuint)));
		}


	}
	void computeSphereTexture(float scale = 1)
	{
		int n = vertices.size() / 3;
		texturecoordinates.resize(2 * n);
		GLfloat x, y, z;
		for (int i = 0; i<n; i++)
		{
			x = vertices[3 * i]; y = vertices[3 * i + 1]; z = vertices[3 * i + 2];
			if (x == 0 && y == 0 && z == 0) continue;
			float l = sqrt(x*x + y*y + z*z);
			x = x / l; y = y / l; z = z / l;

			texturecoordinates[2 * i] = -z >= 0.0f ? atan2(-z, x) / (2 * PI) : (2 * PI + atan2(-z, x)) / (2 * PI);
			texturecoordinates[2 * i + 1] = y >= 0.0f ? acos(y) / PI : (PI - acos(-y)) / PI;

			texturecoordinates[2 * i] *= scale;
			texturecoordinates[2 * i + 1] *= scale;
		}
	}

	void computePlaneTexture()
	{
		int n = vertices.size() / 3;
		texturecoordinates.resize(2 * n);
		int num = sqrt((float)n);

		for (int i = 0; i<num; i++)
			for (int j = 0; j<num; j++){
				texturecoordinates[2 * (i*num + j)] = (float)i / num;
				texturecoordinates[2 * (i*num + j) + 1] = 1.0 - (float)j / num;
			}
	}

	void computeCylinderTexture()
	{
		int n = vertices.size() / 3;
		texturecoordinates.resize(2 * n);
		GLfloat x, y, z;
		for (int i = 0; i<n; i++)
		{
			x = vertices[3 * i]; y = vertices[3 * i + 1]; z = vertices[3 * i + 2];

			texturecoordinates[2 * i] = z;
			if (y >= 0.0f)     texturecoordinates[2 * i + 1] = atan2(y, x) / (PI);
			else if (y<0.0f)  texturecoordinates[2 * i + 1] = (-atan2(y, x)) / (PI);
			//this has problems at the seam, when 1->0 and so interpoltion results in the whole image squeezed between the two border vertices.
			//if ( y>=0.0f )     texturecoordinates[2*i+1] = atan2(  y,  x ) / (2*PI) ;
			//else if ( y<0.0f )  texturecoordinates[2*i+1] = (2*PI + atan2(  y,  x )) / (2*PI) ;
		}
	}

	void computeSphereTextureParts(float scale = 1)
	{
		texturecoordinates.resize(2 * vertices.size() / 3);
		GLfloat x, y, z;
		vector<GLfloat> tmp;
		for (int k = 0; k<parts.size(); k++)
		{
			tmp.clear();
			for (int i = 3 * (parts[k]->start); i<3 * (parts[k]->end); i++)
			{
				int vindex = indices[i];
				tmp.push_back(vertices[3 * vindex]); tmp.push_back(vertices[3 * vindex + 1]); tmp.push_back(vertices[3 * vindex + 2]);
			}
			normalize(tmp);
			for (int i = 3 * (parts[k]->start); i<3 * (parts[k]->end); i++)
			{
				int vindex = indices[i];
				//x = vertices[vindex*3]; y = vertices[vindex*3+1]; z = vertices[vindex*3+2];
				x = tmp[3 * (i - 3 * (parts[k]->start))]; y = tmp[3 * (i - 3 * (parts[k]->start)) + 1]; z = tmp[3 * (i - 3 * (parts[k]->start)) + 2];
				if (x == 0 && y == 0 && z == 0) continue;
				float l = sqrt(x*x + y*y + z*z);
				x = x / l; y = y / l; z = z / l;

				texturecoordinates[2 * vindex] = -z >= 0.0f ? atan2(-z, x) / (2 * PI) : (2 * PI + atan2(-z, x)) / (2 * PI);
				texturecoordinates[2 * vindex + 1] = y >= 0.0f ? acos(y) / PI : (PI - acos(-y)) / PI;

				texturecoordinates[2 * vindex] *= scale;
				texturecoordinates[2 * vindex + 1] *= scale;
			}
			//for (int i=3*parts[k]->start;i<3*parts[k]->end;i++)
			//{
			//	int vindex = indices[i];
			//	x = vertices[vindex*3]; y = vertices[vindex*3+1]; z = vertices[vindex*3+2];
			//	if (x==0&&y==0&&z==0) continue;
			//	float l = sqrt(x*x + y*y + z*z); 
			//	x = x/l; y = y/l; z = z/l;

			//	texturecoordinates[2*vindex] =  -z >= 0.0f ? atan2(  -z,  x ) / (2*PI) : texturecoordinates[2*vindex+1] = (2*PI + atan2(  -z,  x )) / (2*PI) ;
			//	texturecoordinates[2*vindex+1] = y >= 0.0f ? acos(y) / PI : (PI - acos(-y))/ PI;

			//	texturecoordinates[2*vindex] *= scale;
			//	texturecoordinates[2*vindex+1] *= scale;
			//}
		}

	}

	void displayTangents()
	{
		glUseProgram(0);
		int n = vertices.size() / 3;
		double scale = 20.0f;
		glBegin(GL_LINES);
		for (int i = 0; i<n; i++)
		{
			glColor3f(1, 1, 1);
			glVertex3f(vertices[3 * i], vertices[3 * i + 1], vertices[3 * i + 2]);
			glVertex3f(vertices[3 * i] + normals[3 * i] / scale, vertices[3 * i + 1] + normals[3 * i + 1] / scale, vertices[3 * i + 2] + normals[3 * i + 2] / scale);

			glVertex3f(vertices[3 * i], vertices[3 * i + 1], vertices[3 * i + 2]);
			glVertex3f(vertices[3 * i] + tangents[3 * i] / scale, vertices[3 * i + 1] + tangents[3 * i + 1] / scale, vertices[3 * i + 2] + tangents[3 * i + 2] / scale);

			//Computing bitangent by taking cross-product of normal and tangent.
			GLfloat bitangentx = normals[3 * i + 1] * tangents[3 * i + 2] - normals[3 * i + 2] * tangents[3 * i + 1];
			GLfloat bitangenty = normals[3 * i + 2] * tangents[3 * i] - normals[3 * i] * tangents[3 * i + 2];
			GLfloat bitangentz = normals[3 * i] * tangents[3 * i + 1] - normals[3 * i + 1] * tangents[3 * i];

			glVertex3f(vertices[3 * i], vertices[3 * i + 1], vertices[3 * i + 2]);
			glVertex3f(vertices[3 * i] + bitangentx / scale, vertices[3 * i + 1] + bitangenty / scale, vertices[3 * i + 2] + bitangentz / scale);
		}
		glEnd();
		glUseProgram(myshaderprogram);
	}

};
