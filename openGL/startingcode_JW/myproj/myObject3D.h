#pragma once
#include <math.h>
#include <GL/glew.h>
#include <vector>
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
	GLuint buffers[6];

	GLuint myshaderprogram; //the ID of shader program, in case multiple shaders
	myMaterial *material; //contains the color description of the object.
	vector<GLfloat> vertices;
	vector<GLfloat> colors;
	vector<GLuint> indices;
	vector<GLfloat> normals;
	vector<GLfloat> textures;
	vector<GLfloat> tangents;
	

	vector<mySubObject3D *> parts; //contains subparts of the scene.

	//Model matrix
	glm::mat4 model_matrix;
	
	myObject3D(GLuint shader) {
		clear();
		myshaderprogram = shader;
		model_matrix = glm::mat4(1.0f);

		material = (myMaterial*)malloc(sizeof(myMaterial));

	}

	void clear() {
		vertices.clear();
		indices.clear();
		normals.clear();
		textures.clear();
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

	void createObjectBuffers()
	{
		glGenBuffers(4, buffers);

		glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
		glBufferData(GL_ARRAY_BUFFER, vertices.size()*4, &vertices.front(), GL_STATIC_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, buffers[1]);
		glBufferData(GL_ARRAY_BUFFER, indices.size()*4, &indices.front(), GL_STATIC_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, buffers[2]);
		glBufferData(GL_ARRAY_BUFFER, normals.size()*4, &normals.front(), GL_STATIC_DRAW);

		//(GL_ARRAY_BUFFER, buffers[4]);
		//glBufferData(GL_ARRAY_BUFFER, textures.size() * 4, &textures.front(), GL_STATIC_DRAW);

		//glBindBuffer(GL_ARRAY_BUFFER, buffers[5]);
		//glBufferData(GL_ARRAY_BUFFER, tangents.size() * 4, &tangents.front(), GL_STATIC_DRAW);
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
		float du1 = textures[2 * v1] - textures[2 * v0];
		float dv1 = textures[2 * v1 + 1] - textures[2 * v0 + 1];
		float du2 = textures[2 * v2] - textures[2 * v0];
		float dv2 = textures[2 * v2 + 1] - textures[2 * v0 + 1];

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
	void readScene(string filename, bool tonormalize = 0)
	{
		cout << "Start of read scene" << endl;
		clear();
		int i = 0;
		string s, t, u, v, mtlfilename;
		myMaterial current;
		float d;

		ifstream fin(filename);

		parts.push_back(new mySubObject3D(0, 0));

		while (getline(fin, s))
		{
			// << s << endl;
			stringstream myline(s);
			myline >> t;
			if (t == "g")
			{
				parts[i]->end = indices.size() / 3;
				parts[i]->material.setMaterial(current);

				parts.push_back(new mySubObject3D(indices.size() / 3, indices.size() / 3));
				i++;
				myline >> parts[i]->name;
			}
			else if (t == "v")
			{
				myline >> u;
				d = atof((u.substr(0, u.find("/"))).c_str());
				vertices.push_back(d);

				myline >> u;
				d = atof((u.substr(0, u.find("/"))).c_str());
				vertices.push_back(d);

				myline >> u;
				d = atof((u.substr(0, u.find("/"))).c_str());
				vertices.push_back(d);
			}
			else if (t == "mtllib")
			{
				myline >> mtlfilename;
			}
			else if (t == "usemtl")
			{
				myline >> u;
				if (mtlfilename.empty()) { cout << "Error! No material file given.\n"; exit(0); }

				string f = mtlfilename;
				ifstream mtlfin(f);

				while (mtlfin >> v)
				{
					if (v == u)
					{
						while (mtlfin >> v) {
							if (v == "Ns") { mtlfin >> d; current.material_Sh = d; }
							else if (v == "Ka") {
								mtlfin >> d; current.material_Ka[0] = d; mtlfin >> d; current.material_Ka[1] = d;
								mtlfin >> d; current.material_Ka[2] = d;
							}
							else if (v == "Kd") {
								mtlfin >> d; current.material_Kd[0] = d; mtlfin >> d; current.material_Kd[1] = d;
								mtlfin >> d; current.material_Kd[2] = d;
							}
							else if (v == "Ks") {
								mtlfin >> d; current.material_Ks[0] = d; mtlfin >> d; current.material_Ks[1] = d;
								mtlfin >> d; current.material_Ks[2] = d;
							}
							else if (v == "newmtl") break;
						}
						mtlfin.close();
						break;
					}
				}
			}
			else if (t == "s") {}
			else if (t == "f")
			{
				myline >> u;
				int f = atoi((u.substr(0, u.find("/"))).c_str());

				myline >> u;
				int l = atoi((u.substr(0, u.find("/"))).c_str());

				while (myline >> u)
				{
					indices.push_back(f - 1);
					indices.push_back(l - 1);
					l = atoi((u.substr(0, u.find("/"))).c_str());
					indices.push_back(l - 1);
				}
			}
		}
		parts[i]->end = indices.size() / 3;
		parts[i]->material.setMaterial(current);
		if (tonormalize) normalize(vertices);

		cout << "End of read scene" << endl;
	}

	void displayScene(glm::mat4 viewmatrix){
		glUniformMatrix4fv(glGetUniformLocation(myshaderprogram, "mymodel_matrix"), 1, GL_FALSE, &model_matrix[0][0]);
		glm::mat3 normal_matrix = glm::transpose(glm::inverse(glm::mat3(viewmatrix*model_matrix)));
		glUniformMatrix3fv(glGetUniformLocation(myshaderprogram, "mynormal_matrix"), 1, GL_FALSE, &normal_matrix[0][0]);

		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

		glEnableVertexAttribArray(2);
		glBindBuffer(GL_ARRAY_BUFFER, buffers[2]);
		glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, 0);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffers[1]);
		glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);

		glUniformMatrix4fv(glGetUniformLocation(myshaderprogram, "mymodel_matrix"),
			1, GL_FALSE, &model_matrix[0][0]);

		for (int i = 0; i<parts.size(); i++) {
			parts[i]->material.activateMaterial();
			/*parts[i]->texture.activateTexture(GL_TEXTURE1);
			parts[i]->bump.activateTexture(GL_TEXTURE2);
			parts[i]->cubemap.activateTexture(GL_TEXTURE3, GL_TEXTURE_CUBE_MAP);
			*/
			glUniform4fv(glGetUniformLocation(myshaderprogram, "mymodel_kd"), 1, parts[i]->material.material_Kd);
			glUniform4fv(glGetUniformLocation(myshaderprogram, "mymodel_ka"), 1, parts[i]->material.material_Ka);
			glUniform4fv(glGetUniformLocation(myshaderprogram, "mymodel_ks"), 1, parts[i]->material.material_Ks);
			glUniform1f(glGetUniformLocation(myshaderprogram, "mymodel_ns"), parts[i]->material.material_Sh);
			glDrawElements(GL_TRIANGLES, (parts[i]->end - parts[i]->start) * 3, GL_UNSIGNED_INT, (GLvoid*)(sizeof(GLuint) * parts[i]->start * 3));
		}


	}
};
