#include <math.h>
#include <GL/glew.h>
#include <vector>
#include <string>
#include <fstream>
#include "vector3d.h"
#include "myMaterial.h"

#define PI 3.14159265

using namespace std;

class myObject3D
{
public:
	GLuint buffers[6];

	GLuint myshaderprogram; //the ID of shader program, in case multiple shaders
	myMaterial *material; //contains the color description of the object.
	vector<GLfloat> vertices; //vertices array.
	vector<GLuint> indices; //indices array.
	vector<GLfloat> normals; //normals array.

	//Model matrix
	glm::mat4 model_matrix;
	
	myObject3D() {
		model_matrix = glm::mat4(1.0f);

		material = (myMaterial*)malloc(sizeof(myMaterial));

	}

	void clear() {
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
		glGenBuffers(3, buffers);

		glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
		glBufferData(GL_ARRAY_BUFFER, vertices.size()*4, &vertices.front(), GL_STATIC_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, buffers[1]);
		glBufferData(GL_ARRAY_BUFFER, indices.size()*4, &indices.front(), GL_STATIC_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, buffers[2]);
		glBufferData(GL_ARRAY_BUFFER, normals.size()*4, &normals.front(), GL_STATIC_DRAW);
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
};
