#include <math.h>
#include <GL/glew.h>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include "vector3d.h"
#include "myMaterial.h"
#include "mySubObject3D.h"
#define PI 3.14159265


using namespace std;
class myObject3D
{
public:
  GLuint buffers[6]; //object buffers to store vertices, indices, normals etc.
  GLuint myshaderprogram; //the ID of shader program, in case multiple shaders
  myMaterial *material; //contains the color description of the object.
  vector<GLfloat> vertices; //vertices array.
  vector<GLuint> indices; //indices array.
  vector<GLfloat> normals; //normals array.
  vector<GLfloat> colors;
  vector<GLfloat> textures;
  vector<GLfloat> tangents;


  vector<mySubObject3D *> parts; //contains subparts of the scene.

  //Model matrix
  glm::mat4 model_matrix;

  myObject3D() {
    model_matrix = glm::mat4(1.0f);

    material = (myMaterial*)malloc(sizeof(myMaterial));
    material->material_Ks[0]= 1;material->material_Ks[1]=1;material->material_Ks[2] =1; material->material_Ks[3] = 0;
    material->material_Kd[0]= 1;material->material_Kd[1]=1;material->material_Kd[2] =1; material->material_Kd[3] = 0;


  }

  void clear() {
	  	vertices.clear();
		colors.clear();
		indices.clear();
		normals.clear();
		textures.clear();
		parts.clear();
		tangents.clear();
  }

  void setShader(GLuint & s)
	{
		myshaderprogram = s;
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
	//ok
  }

  

  void readScene2(string filename, bool tonormalize = 0)
  {	
	  //lol
  }
	void readScene(string filename, bool tonormalize = 0)
	{
		clear();
		int i = 0;
		string s, t, u, v, mtlfilename;
		myMaterial current;
		float d;

		ifstream fin(filename);

		parts.push_back( new mySubObject3D(0,0) );

		while ( getline(fin, s) )
		{
			//cout << s << endl;
			stringstream myline(s);
			myline >> t;
			if ( t == "g" )
			{
				parts[i]->end = indices.size()/3;
				parts[i]->material.setMaterial(current);

				parts.push_back( new mySubObject3D(indices.size()/3, indices.size()/3) );
				i++;
				myline >> parts[i]->name;
			}
			else if (t == "v")
			{
				myline >> u;
				d = atof( (u.substr(0, u.find("/"))).c_str());
				vertices.push_back(d);

				myline >> u;
				d = atof( (u.substr(0, u.find("/"))).c_str());
				vertices.push_back(d);

				myline >> u;
				d = atof( (u.substr(0, u.find("/"))).c_str());
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

				while ( mtlfin >> v )
				{
					if (v == u)
					{
						while (mtlfin >> v) {
							if (v == "Ns") { mtlfin >> d; current.material_Sh = d; }
							else if (v == "Ka") { mtlfin >> d; current.material_Ka[0] = d; mtlfin >> d; current.material_Ka[1] = d; 
							mtlfin >> d; current.material_Ka[2] = d;  }
							else if (v == "Kd") { mtlfin >> d; current.material_Kd[0] = d; mtlfin >> d; current.material_Kd[1] = d; 
							mtlfin >> d; current.material_Kd[2] = d;   }
							else if (v == "Ks") { mtlfin >> d; current.material_Ks[0] = d; mtlfin >> d; current.material_Ks[1] = d; 
							mtlfin >> d; current.material_Ks[2] = d;  }
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
				int f = atoi( (u.substr(0, u.find("/"))).c_str());

				myline >> u;
				int l = atoi( (u.substr(0, u.find("/"))).c_str());

				while ( myline >> u )
				{
					indices.push_back( f-1 );
					indices.push_back( l-1 );
					l = atoi( (u.substr(0, u.find("/"))).c_str());
					indices.push_back( l-1 );
				}
			}
		}
		parts[i]->end = indices.size()/3;
		parts[i]->material.setMaterial(current);
		if (tonormalize) normalize(vertices);
	}
	
	void computeTangent(int v0, int v1, int v2, float & x, float & y, float & z)
	{
		float du1 = textures[2*v1] - textures[2*v0] ;
		float dv1 = textures[2*v1+1] - textures[2*v0+1] ;
		float du2 = textures[2*v2] - textures[2*v0] ;
		float dv2 = textures[2*v2+1] - textures[2*v0+1];

		float f = 1.0f / (du1 * dv2 - du2 * dv1);
		if ( (du1*dv2 - du2*dv1) == 0){
			x = y = z = 0; return;
		}

		float e1x = vertices[3*v1] - vertices[3*v0];
		float e1y = vertices[3*v1+1] - vertices[3*v0+1];
		float e1z = vertices[3*v1+2] - vertices[3*v0+2];

		float e2x = vertices[3*v2] - vertices[3*v0];
		float e2y = vertices[3*v2+1] - vertices[3*v0+1];
		float e2z = vertices[3*v2+2] - vertices[3*v0+2];

		x = f * ( dv2 * e1x - dv1 * e2x );
		y = f * ( dv2 * e1y - dv1 * e2y );
		z = f * ( dv2 * e1z - dv1 * e2z );
	}
	void computeTangents( )
	{
		int i, j, k;
		GLfloat x1, y1, z1;

		int n = vertices.size()/3;
		int m = indices.size()/3;

		tangents.resize(3*n);
		int *incidences = new int[n];
		for (i=0;i<3*n;i++) tangents[i] = 0.0;
		for (i=0;i<n;i++) incidences[i] = 0;

		for (j=0;j<m;j++)
		{
			computeTangent(indices[3*j], indices[3*j+1], indices[3*j+2], x1, y1, z1);
			tangents[3*indices[3*j]] += x1; tangents[3*indices[3*j]+1] += y1; tangents[3*indices[3*j]+2] += z1;
			tangents[3*indices[3*j+1]] += x1; tangents[3*indices[3*j+1]+1] += y1; tangents[3*indices[3*j+1]+2] += z1;
			tangents[3*indices[3*j+2]] += x1; tangents[3*indices[3*j+2]+1] += y1; tangents[3*indices[3*j+2]+2] += z1;
			incidences[indices[3*j]]++; incidences[indices[3*j+1]]++; incidences[indices[3*j+2]]++;
		}
		for (i=0;i<n;i++) {
			float l = sqrt( tangents[3*i]*tangents[3*i] + tangents[3*i+1]*tangents[3*i+1] + tangents[3*i+2]*tangents[3*i+2] );
			tangents[3*i] /= l; tangents[3*i+1] /= l; tangents[3*i+2] /= l;
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
		glGenBuffers(6, buffers);

		glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
		glBufferData(GL_ARRAY_BUFFER, vertices.size()*4, &vertices.front(), GL_STATIC_DRAW);

		/*glBindBuffer(GL_ARRAY_BUFFER, buffers[1]) ; 
		if (!colors.empty()) glBufferData(GL_ARRAY_BUFFER, n*3*4, &colors.front(), GL_STATIC_DRAW);*/

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffers[2]);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size()*4, &indices.front(), GL_STATIC_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, buffers[3]);
		glBufferData(GL_ARRAY_BUFFER, normals.size()*4, &normals.front(), GL_STATIC_DRAW);

	    glBindBuffer(GL_ARRAY_BUFFER, buffers[4]) ; 
		glBufferData(GL_ARRAY_BUFFER, textures.size()*4, &textures.front(), GL_STATIC_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, buffers[5]) ; 
		glBufferData(GL_ARRAY_BUFFER, tangents.size()*4, &tangents.front(), GL_STATIC_DRAW);
  }

  void displayObject(int c = -1)
  {
	  	GLuint tangent_loc = glGetAttribLocation(myshaderprogram,"tangent");
	
	
		glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
		glVertexPointer(3, GL_FLOAT, 0, 0);
		glEnableClientState(GL_VERTEX_ARRAY);

		//glBindBuffer(GL_ARRAY_BUFFER, buffers[1]);
		//glColorPointer(3, GL_FLOAT, 0, 0) ; 
		//glEnableClientState(GL_COLOR_ARRAY) ; 

		glBindBuffer(GL_ARRAY_BUFFER, buffers[3]);
    	glNormalPointer(GL_FLOAT, 0, 0) ; 
		glEnableClientState(GL_NORMAL_ARRAY) ;

		glBindBuffer(GL_ARRAY_BUFFER, buffers[4]);
		glTexCoordPointer(2,GL_FLOAT,0,0) ; 
		glEnableClientState(GL_TEXTURE_COORD_ARRAY) ; 

		glBindBuffer(GL_ARRAY_BUFFER, buffers[5]);
		glVertexAttribPointer(tangent_loc, 3, GL_FLOAT, 0, 0, 0);
		glEnableVertexAttribArray(tangent_loc);
		
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffers[2]);

		if (c != -1) {
			c = c % parts.size();
			parts[c]->material.activateMaterial();
			parts[c]->texture.activateTexture(GL_TEXTURE1);
			parts[c]->bump.activateTexture(GL_TEXTURE2);
			parts[c]->cubemap.activateTexture(GL_TEXTURE3, GL_TEXTURE_CUBE_MAP);
			glDrawElements(GL_TRIANGLES, (parts[c]->end-parts[c]->start)*3, GL_UNSIGNED_INT, (GLvoid*)(sizeof(GLuint) * parts[c]->start * 3)) ; 	
			return;
		}
		for (int i=0;i<parts.size();i++) {
			parts[i]->material.activateMaterial();
			parts[i]->texture.activateTexture(GL_TEXTURE1);
			parts[i]->bump.activateTexture(GL_TEXTURE2);
			parts[i]->cubemap.activateTexture(GL_TEXTURE3, GL_TEXTURE_CUBE_MAP);
			glDrawElements(GL_TRIANGLES, (parts[i]->end-parts[i]->start)*3, GL_UNSIGNED_INT, (GLvoid*)(sizeof(GLuint) * parts[i]->start * 3)) ; 
		}
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

  void normalize(vector<GLfloat> & myvertices)
	{
		if (myvertices.size() < 4) return;
		int i;
		int tmpxmin = 0, tmpymin = 0, tmpzmin = 0, tmpxmax = 0, tmpymax = 0, tmpzmax = 0;
		
		int n = myvertices.size()/3;
		
		for (i=0;i<n;i++) {
			if (myvertices[3*i] < myvertices[3*tmpxmin]) tmpxmin = i;
			if (myvertices[3*i] > myvertices[3*tmpxmax]) tmpxmax = i;

			if (myvertices[3*i+1] < myvertices[3*tmpymin+1]) tmpymin = i;
			if (myvertices[3*i+1] > myvertices[3*tmpymax+1]) tmpymax = i;

			if (myvertices[3*i+2] < myvertices[3*tmpzmin+2]) tmpzmin = i;
			if (myvertices[3*i+2] > myvertices[3*tmpzmax+2]) tmpzmax = i;
		}

		double xmin = myvertices[3*tmpxmin], xmax = myvertices[3*tmpxmax], 
			   ymin = myvertices[3*tmpymin+1], ymax = myvertices[3*tmpymax+1], 
			   zmin = myvertices[3*tmpzmin+2], zmax = myvertices[3*tmpzmax+2];

		double scale = (xmax-xmin) <= (ymax-ymin) ? (xmax-xmin) : (ymax-ymin);
		scale = scale >= (zmax-zmin) ? scale : (zmax-zmin);

		for (i=0;i<n;i++) {
			myvertices[3*i] -= (xmax+xmin)/2;
			myvertices[3*i+1] -= (ymax+ymin)/2;
			myvertices[3*i+2] -= (zmax+zmin)/2;

			myvertices[3*i] /= scale;
			myvertices[3*i+1] /= scale;
			myvertices[3*i+2] /= scale;
		}
	}
  
	void setMaterial(GLfloat ka_r, GLfloat ka_g, GLfloat ka_b, GLfloat ka_a, GLfloat kd_r, GLfloat kd_g, GLfloat kd_b, GLfloat kd_a, GLfloat ks_r, GLfloat ks_g, GLfloat ks_b, GLfloat ks_a, GLfloat s) 
	{
		for (int i=0;i<parts.size();i++)
			parts[i]->material.setMaterial(ka_r,  ka_g,  ka_b,  ka_a,  kd_r,  kd_g,  kd_b,  kd_a, ks_r,  ks_g,  ks_b,  ks_a,  s);
	}

	void setTexture(myTexture *t)
	{
		for (int i=0;i<parts.size();i++)
			parts[i]->texture = *t;
	}

	void setBump(myTexture *t)
	{
		for (int i=0;i<parts.size();i++)
			parts[i]->bump = *t;
	}

	void setCubemap(myTexture *t)
	{
		for (int i=0;i<parts.size();i++)
			parts[i]->cubemap = *t;
	}
	
};
