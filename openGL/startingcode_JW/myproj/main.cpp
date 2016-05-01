#include "stdafx.h"
#include <fstream>
#include <GL/glew.h>
#include <freeglut.h>
#include <string>
#include <sstream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp> 
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <cstdlib>
#include <iostream>
#include <ctime>

using namespace std;

#include "shaders.h"

#include "point3d.h"
#include "vector3d.h"
#include "myObject3D.h"
#include "myLight.h"
#include <math.h>



// width and height of the window.
int Glut_w = 600, Glut_h = 400; 

//Variables and their values for the camera setup.
myPoint3D camera_eye(0,3,-1);
myVector3D camera_up(0,1,0);
myVector3D camera_forward (0,0,1);
myVector3D camera_side(1, 0, 0);

float fovy = 90;
float zNear = 0.2;
float zFar = 6000;

int button_pressed = 0; // 1 if a button is currently being pressed.
int GLUTmouse[2] = { 0, 0 };

GLuint vertexshader, fragmentshader, shaderprogram1; // shaders

GLuint renderStyle = 0;			GLuint renderStyle_loc;
GLuint projection_matrix_loc;
GLuint view_matrix_loc;
GLuint normal_matrix_loc;
GLuint buffers[6];

vector<GLfloat> vertices;
vector<GLfloat> normals;
vector<GLuint> indices;

myObject3D *obj1;
myObject3D *obj2;
myObject3D *obj3;
myObject3D *obj4;

myLight *lights;
int nbLight;


//This function is called when a mouse button is pressed.
void mouse(int button, int state, int x, int y)
{
  // Remember button state 
  button_pressed = (state == GLUT_DOWN) ? 1 : 0;

   // Remember mouse position 
  GLUTmouse[0] = x;
  GLUTmouse[1] = Glut_h - y;
}

//This function is called when the mouse is dragged.
void mousedrag(int x, int y)
{
  // Invert y coordinate
  y = Glut_h - y;

  //change in the mouse position since last time
  int dx = x - GLUTmouse[0];
  int dy = y - GLUTmouse[1];

  GLUTmouse[0] = x;
  GLUTmouse[1] = y;

  if (dx == 0 && dy == 0) return;
  if (button_pressed == 0) return;

  double vx = (double) dx / (double) Glut_w;
  double vy = (double) dy / (double) Glut_h;
  double theta = 4.0 * (fabs(vx) + fabs(vy));

  myVector3D camera_right = camera_forward.crossproduct(camera_up);
  camera_right.normalize();

  myVector3D tomovein_direction = -camera_right*vx + -camera_up*vy;

  myVector3D rotation_axis = tomovein_direction.crossproduct(camera_forward);
  rotation_axis.normalize();
  
  camera_forward.rotate(rotation_axis, theta);
  
  camera_up.rotate(rotation_axis, theta);
  camera_eye.rotate(rotation_axis, theta);
 
  camera_up.normalize();
  camera_forward.normalize();

  glutPostRedisplay();
}

void mouseWheel(int button, int dir, int x, int y)
{
	if (dir > 0)
		camera_eye += camera_forward * 0.02;
	else
		camera_eye += -camera_forward * 0.02;
	glutPostRedisplay();
}

//This function is called when a key is pressed.
void keyboard(unsigned char key, int x, int y) {
	switch(key) {
	case 27:  // Escape to quit
		exit(0) ;
        break ;
	case 'f':
		renderStyle = (renderStyle+1)%2;
		glUniform1i(renderStyle_loc, renderStyle) ; 
		break;
	case 'r':
		camera_eye = myPoint3D(0,0,2);
		camera_up = myVector3D(0,1,0);
		camera_forward = myVector3D(0,0,-1);
		break;
	case  'z': 
		camera_eye.X -= 0.1 * sin(180 * PI / 180.0);
		camera_eye.Z -= 0.1 * cos(180 * PI / 180.0);
		obj4->translate(-0.1 * sin(180 * PI / 180.0), 0, -0.1 * cos(180 * PI / 180.0));
		break;
	case  's':
		camera_eye.X += 0.1 * sin(180 * PI / 180.0);
		camera_eye.Z += 0.1 * cos(180 * PI / 180.0);
		obj4->translate(0.1 * sin(180 * PI / 180.0), 0, 0.1 * cos(180 * PI / 180.0));
		break;
	case  'q':
		camera_eye += myVector3D(0.1, 0, 0);
		obj4->translate(0.1, 0, 0);
		break;
	case  'd':
		camera_eye += myVector3D(-0.1, 0, 0);
		obj4->translate(-0.1, 0, 0);
		break;
	}
	glutPostRedisplay();
}

//This function is called when an arrow key is pressed.
void keyboard2(int key, int x, int y) {
	switch(key) {
	case GLUT_KEY_UP:
		//camera_eye += camera_forward*0.1;
		camera_side.normalize();
		camera_forward.rotate(camera_side, 0.1);
		camera_forward.normalize();
		break;
	case GLUT_KEY_DOWN:
		//camera_eye += -camera_forward*0.1;
		camera_side.normalize();
		camera_forward.rotate(camera_side, -0.1);
		camera_forward.normalize();
		break;
	case GLUT_KEY_LEFT:
		camera_up.normalize();
		camera_forward.rotate(camera_up, 0.1);
		camera_forward.normalize();
		break;
	case GLUT_KEY_RIGHT:
		camera_up.normalize();
		camera_forward.rotate(camera_up, -0.1);
		camera_forward.normalize();
		break;
	}
	glutPostRedisplay();
}

void reshape(int width, int height){
	Glut_w = width;
	Glut_h = height;
	glm::mat4 projection_matrix = glm::perspective(fovy, Glut_w/(float)Glut_h, zNear, zFar);
	glUniformMatrix4fv(projection_matrix_loc, 1, GL_FALSE, &projection_matrix[0][0]);
	glViewport(0, 0, Glut_w, Glut_h);
}

//This function is called to display objects on screen.
void display() 
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0, Glut_w, Glut_h);

	glm::mat4 projection_matrix = 
		glm::perspective(fovy, Glut_w/(float)Glut_h, zNear, zFar);
	glUniformMatrix4fv(projection_matrix_loc, 1, GL_FALSE, &projection_matrix[0][0]);

	glm::mat4 view_matrix = 
		glm::lookAt(glm::vec3(camera_eye.X, camera_eye.Y, camera_eye.Z), 
					glm::vec3(camera_eye.X + camera_forward.dX, camera_eye.Y + camera_forward.dY, camera_eye.Z + camera_forward.dZ), 
					glm::vec3(camera_up.dX, camera_up.dY, camera_up.dZ));
	glUniformMatrix4fv(view_matrix_loc, 1, GL_FALSE, &view_matrix[0][0]);

	glm::mat3 normal_matrix = glm::transpose(glm::inverse(glm::mat3(view_matrix)));
	glUniformMatrix4fv(normal_matrix_loc, 1, GL_FALSE, &normal_matrix[0][0]);




	
	// 3.6. Display the light
	// - point for point light
	// - segment for direction light
	glUniform1i(renderStyle_loc, 2);
	for(int i = 0; i<nbLight; i++){
		lights[i].drawLight(shaderprogram1);
	}
	
	//obj1->displayObject(shaderprogram1,view_matrix);
	//obj2->displayObject(shaderprogram1,view_matrix);
	obj3->displayScene(view_matrix);
	obj4->displayScene(view_matrix);

	glFlush();
}

//This function is called from the main to initalize everything.
void init()
{
    vertexshader = initshaders(GL_VERTEX_SHADER, "shaders/light.vert.glsl") ;
    fragmentshader = initshaders(GL_FRAGMENT_SHADER, "shaders/light.frag.glsl") ;
    shaderprogram1 = initprogram(vertexshader, fragmentshader);

	renderStyle_loc = glGetUniformLocation(shaderprogram1, "myrenderStyle") ;
	glUniform1i(renderStyle_loc, renderStyle);

	projection_matrix_loc = glGetUniformLocation(shaderprogram1, "myprojection_matrix");
	view_matrix_loc = glGetUniformLocation(shaderprogram1, "myview_matrix");
	normal_matrix_loc = glGetUniformLocation(shaderprogram1, "mynormal_matrix");





	/*obj1 = new myObject3D();
	obj1->readMesh("apple.obj");
	obj1->computeNormals();
	obj1->createObjectBuffers();*/

	// Second object
	/*obj2 = new myObject3D();
	obj2->readMesh("hand.obj");
	obj2->computeNormals();
	obj2->createObjectBuffers();
	obj2->translate(4,0,0);*/
	

	
	// Maya Scene
	obj3 = new myObject3D(shaderprogram1);
	//obj3->readScene("test.obj");
	//obj3->readScene("TheCarnival.obj");
	obj3->readScene("portal2.obj");
	if (obj3->normals.size() == 0)
		obj3->computeNormals();
	obj3->createObjectBuffers();


	// Personnage iron man
	obj4 = new myObject3D(shaderprogram1);
	obj4->readScene("IronMan2.obj");
	if (obj4->normals.size() == 0)
		obj4->computeNormals();
	obj4->createObjectBuffers();
	obj4->scale(0.008, 0.008, 0.008);
	obj4->translate(0, 1, -1);
	//obj4->rotate(0,1,0,180);

	// 3.5
	nbLight = 1;
	int l1 = 0, l2 = 0, l3 = 2;

	lights = (myLight*)malloc(nbLight*sizeof(myLight));

	lights[l1].color[0] = 1; lights[l1].color[1] = 1; lights[l1].color[2] = 1; lights[l1].color[3] = 1;
	lights[l1].position[0] = 0; lights[l1].position[1] = 4; lights[l1].position[2] = 0; lights[l1].position[3] = 0;
	lights[l1].direction[0] = 0; lights[l1].direction[1] = -3; lights[l1].direction[2] = 0; lights[l1].direction[3] = 0;
	lights[l1].type = 0; // point light
	/*
	lights[l2].color[0] = 1; lights[l2].color[1] = 1; lights[l2].color[2] = 1; lights[l2].color[3] = 1;
	lights[l2].position[0] = 8; lights[l2].position[1] = 20; lights[l2].position[2] = -10; lights[l2].position[3] = 0;
	lights[l2].direction[0] = 0; lights[l2].direction[1] = -1; lights[l2].direction[2] = 0; lights[l2].direction[3] = 0;
	lights[l2].type = 1; // direction light
	

	lights[l3].color[0] = 1; lights[l3].color[1] = 1; lights[l3].color[2] = 1; lights[l3].color[3] = 0; 
	lights[l3].position[0] = -1; lights[l3].position[1] = 2; lights[l3].position[2] = 0; lights[l3].position[3] = 0; 
	lights[l3].direction[0] = 0; lights[l3].direction[1] = 2; lights[l3].direction[2] = 0; lights[l3].direction[3] = 0; 
	lights[l3].type = 2; // spotlight
	*/
	glUniform1i(glGetUniformLocation(shaderprogram1, "nbLights"), nbLight);

	// 3.5

	vector<GLfloat> colors;
	vector<GLfloat> direction;
	vector<GLfloat> position;
	vector<GLint> type;

	for (int i = 0; i<nbLight; i++){
		colors.insert(colors.end(), lights[i].color, lights[i].color + 4);
		direction.insert(direction.end(), lights[i].direction, lights[i].direction + 4);
		position.insert(position.end(), lights[i].position, lights[i].position + 4);
		type.push_back(lights[i].type);
	}


	glUniform4fv(glGetUniformLocation(shaderprogram1, "light_colors"), nbLight, &(colors.front()));
	glUniform4fv(glGetUniformLocation(shaderprogram1, "light_position"), nbLight, &(position.front()));
	glUniform3fv(glGetUniformLocation(shaderprogram1, "light_direction"), nbLight, &(direction.front()));
	glUniform1iv(glGetUniformLocation(shaderprogram1, "light_type"), nbLight, &(type.front()));


	//glUniform1i(glGetUniformLocation(shaderprogram1, "tex"), 8);


	//glDisable(GL_TEXTURE_2D);
	glClearColor(0.4, 0.4, 0.4, 0);

}


int main(int argc, char* argv[]) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA | GLUT_DEPTH );
	glutCreateWindow("My OpenGL Application");
	   
	glewInit() ; 
	glutReshapeWindow(Glut_w, Glut_h);
	
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(keyboard2);
	glutMotionFunc(mousedrag) ;
	glutMouseFunc(mouse) ;
	glutMouseWheelFunc(mouseWheel);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS) ;
	
	init();

	glutMainLoop();
	return 0;
}