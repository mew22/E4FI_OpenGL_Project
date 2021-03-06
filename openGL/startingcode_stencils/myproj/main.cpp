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

//Variables for window settings
int Glut_w = 600, Glut_h = 400; 
float fovy = 90;
float zNear = 0.2;
float zFar = 6000;

//Variables for mouse 
int button_pressed = 0; // 1 if a button is currently being pressed.
int GLUTmouse[2] = { 0, 0 };

//Variable for shader
GLuint vertexshader, fragmentshader, shaderprogram1;


//Variable for shader location
GLuint renderStyle = 0;			
GLuint renderStyle_loc;
GLuint projection_matrix_loc;
GLuint view_matrix_loc;
GLuint normal_matrix_loc;
GLuint buffers[6];

//Variables for objects models
myObject3D *obj1;
myObject3D *obj2;
myObject3D *level;
myObject3D *character;


//Variables for lighting
myLight *lights;
int nbLight;

//Variables and their values for the camera setup.
myPoint3D camera_eye(0, 4, -1);
myVector3D camera_up(0, 1, 0);
myVector3D camera_forward(0, 0, 1);
myVector3D camera_side(1, 0, 0);
myVector3D camera_right;
myVector3D camera_target;
float sensitivity = 0.2;
float speedX = 0.1, speedY = 0, speedZ = 0.1;
float teta = 180, phi = 0;
float r_temp = 1;
float arm_angle = 0;

//Variables for portal creation
bool portalActivated[2] = { false, false };
myObject3D portals[2];
myVector3D portalsDirection[2];
int portalsDirectionIteraction[2] = { 0,0 };
int portalDistance = 10;


//This function is called when a mouse button is pressed.
void mouse(int button, int state, int x, int y)
{
  // Remember button state 
  button_pressed = (state == GLUT_DOWN) ? 1 : 0;

   /*
  GLUTmouse[0] = x;
  GLUTmouse[1] = Glut_h - y;
  */

  if (button_pressed == 1 && button == GLUT_LEFT_BUTTON) {
	  cout << "portal1 created" << endl;

	  portals[0].model_matrix = character->model_matrix;
	  portalsDirection[0] = camera_forward;
	  portalsDirectionIteraction[0] = portalDistance;
	  portalActivated[0] = true;

  }
  if (button_pressed == 1 && button == GLUT_RIGHT_BUTTON) {
	  cout << "portal2 created" << endl;

	  portals[1].model_matrix = character->model_matrix;
	  portalsDirection[1] = camera_forward;
	  portalsDirectionIteraction[1] = portalDistance;
	  portalActivated[1] = true;
  }
}

//This function is called when the mouse is dragged.
void mousedrag(int x, int y)
{
  // Invert y coordinate
  //y = Glut_h - y;

  //change in the mouse position since last time


  int dx = x - GLUTmouse[0];
  int dy = y - GLUTmouse[1];

  GLUTmouse[0] = x;
  GLUTmouse[1] = y;



  teta -= dx * sensitivity;
  
  phi -= dy * sensitivity;
  if (phi > 89)
	  phi = 89;
  else if (phi < -89)
	  phi = -89;


  if (dx == 0 && dy == 0) return;
  /*if (button_pressed == 0) return;*/

  //float rayon = sqrt(pow(dx, 2) + pow(dy, 2));
  r_temp = cos(phi * PI / 180);
  camera_forward.dX = sin(teta * PI / 180.0) * r_temp;
  camera_forward.dY = sin(phi * PI / 180.0);
  camera_forward.dZ = cos(teta * PI / 180.0) * r_temp;

  if (character != nullptr) {
	  character->translate(-camera_eye.X, 0, -camera_eye.Z);
	  character->rotate(0, 1, 0, -arm_angle);
	  character->rotate(0, 1, 0, teta);
	  character->translate(camera_eye.X, 0, camera_eye.Z);
	  arm_angle = teta;
  }
 

  double vx = (double) dx / (double) Glut_w;
  double vy = (double) dy / (double) Glut_h;
  //double theta = 4.0 * (fabs(vx) + fabs(vy));

  camera_right = camera_forward.crossproduct(camera_up);
  camera_right.normalize();

  camera_target.dX = camera_eye.X + camera_forward.dX;
  camera_target.dY = camera_eye.Y + camera_forward.dY;
  camera_target.dZ = camera_eye.Z + camera_forward.dZ;

  /*myVector3D tomovein_direction = -camera_right*vx + -camera_up*vy;

  myVector3D rotation_axis = tomovein_direction.crossproduct(camera_forward);
  rotation_axis.normalize();
  
  camera_forward.rotate(rotation_axis, theta);
  
  camera_up.rotate(rotation_axis, theta);
  camera_eye.rotate(rotation_axis, theta);
 
  camera_up.normalize();
  camera_forward.normalize();*/

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
		camera_eye.X += speedX * camera_forward.dX;
		camera_eye.Z += speedZ * camera_forward.dZ;

		character->translate(speedX * camera_forward.dX, 0, speedZ * camera_forward.dZ);
		break;
	case  's':
		camera_eye.X -= speedX * camera_forward.dX;
		camera_eye.Z -= speedZ * camera_forward.dZ;

		character->translate(-speedX * camera_forward.dX, 0, -speedZ * camera_forward.dZ);
		break;
	case  'q':
		//camera_eye.X -= speedX * cos(teta * PI / 180.0);
		//camera_eye.Z += speedZ * sin(teta * PI / 180.0);

		camera_eye.X -= speedX * camera_right.dX;
		camera_eye.Z -= speedZ * camera_right.dZ;

		character->translate(-speedX * camera_right.dX, 0, -speedZ * camera_right.dZ);
		break;
	case  'd':
		//camera_eye.X += speedX * cos(teta * PI / 180.0);
		//camera_eye.Z -= speedZ * sin(teta * PI / 180.0);

		camera_eye.X += speedX * camera_right.dX;
		camera_eye.Z += speedZ * camera_right.dZ;

		character->translate(speedX * camera_right.dX, 0, speedZ * camera_right.dZ);
		break;
	}
	camera_target.dX = camera_eye.X + camera_forward.dX;
	camera_target.dY = camera_eye.Y + camera_forward.dY;
	camera_target.dZ = camera_eye.Z + camera_forward.dZ;
	glutPostRedisplay();
}

//This function is called when an arrow key is pressed.
void keyboard2(int key, int x, int y) {
	//perso_cam.deplacer(souris, GLUTmouse[0], GLUTmouse[1], key);
	switch(key) {
	case GLUT_KEY_UP:
		//camera_eye += camera_forward*0.1;
		/*camera_side.normalize();
		camera_forward.rotate(camera_side, 0.1);
		camera_forward.normalize();*/

		phi += 5;
		if (phi > 89)
			phi = 89;
		else if (phi < -89)
			phi = -89;
		camera_forward.dY = sin(phi * PI / 180.0);

		break;
	case GLUT_KEY_DOWN:
		//camera_eye += -camera_forward*0.1;
		/*camera_side.normalize();
		camera_forward.rotate(camera_side, -0.1);
		camera_forward.normalize();*/

		phi -= 5;
		if (phi > 89)
			phi = 89;
		else if (phi < -89)
			phi = -89;
		camera_forward.dY = sin(phi * PI / 180.0);
		

		break;
	case GLUT_KEY_LEFT:
		/*camera_up.normalize();
		camera_forward.rotate(camera_up, 0.1);
		//character->rotate(0, 1, 0, 0.1 * 180 / PI);
		camera_forward.normalize();*/

		r_temp = cos(phi * PI / 180);
		teta += 5;
		camera_forward.dX = sin(teta * PI / 180.0) * r_temp;
		camera_forward.dZ = cos(teta * PI / 180.0) * r_temp;

		if (character != nullptr) {
			character->translate(-camera_eye.X, 0, -camera_eye.Z);
			character->rotate(0, 1, 0, -arm_angle);
			character->rotate(0, 1, 0, teta);
			character->translate(camera_eye.X, 0, camera_eye.Z);
			arm_angle = teta;
		}
		break;

	case GLUT_KEY_RIGHT:
		/*camera_up.normalize();
		camera_forward.rotate(camera_up, -0.1);
		//character->rotate(0, 1, 0, -0.1 * 180 / PI);
		camera_forward.normalize();*/

		r_temp = cos(phi * PI / 180);
		teta -= 5;
		camera_forward.dX = sin(teta * PI / 180.0) * r_temp;
		camera_forward.dZ = cos(teta * PI / 180.0) * r_temp;

		if (character != nullptr) {
			character->translate(-camera_eye.X, 0, -camera_eye.Z);
			character->rotate(0, 1, 0, -arm_angle);
			character->rotate(0, 1, 0, teta);
			character->translate(camera_eye.X, 0, camera_eye.Z);
			arm_angle = teta;
		}
		break;
	}
	camera_right = camera_forward.crossproduct(camera_up);
	camera_right.normalize();

	camera_target.dX = camera_eye.X + camera_forward.dX;
	camera_target.dY = camera_eye.Y + camera_forward.dY;
	camera_target.dZ = camera_eye.Z + camera_forward.dZ;
	glutPostRedisplay();
}

void reshape(int width, int height){
	Glut_w = width;
	Glut_h = height;
	glm::mat4 projection_matrix = glm::perspective(fovy, Glut_w/(float)Glut_h, zNear, zFar);
	glUniformMatrix4fv(projection_matrix_loc, 1, GL_FALSE, &projection_matrix[0][0]);
	glViewport(0, 0, Glut_w, Glut_h);
}

glm::mat4 portal_view(glm::mat4 orig_view, myObject3D& src, myObject3D& dst) {
	glm::mat4 mv = orig_view * src.model_matrix;
	glm::mat4 portal_cam =
		// 3. transformation from source portal to the camera - it's the
		//    first portal's ModelView matrix:
		mv
		// 2. object is front-facing, the camera is facing the other way:
		* glm::rotate(glm::mat4(1.0), 180.0f, glm::vec3(0.0, 1.0, 0.0))
		// 1. go the destination portal; using inverse, because camera
		//    transformations are reversed compared to object
		//    transformations:
		* glm::inverse(dst.model_matrix)
		;
	return portal_cam;
}

//This function is called to display objects on screen.
void display() 
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	glViewport(0, 0, Glut_w, Glut_h);

	glm::mat4 projection_matrix = 
		glm::perspective(fovy, Glut_w/(float)Glut_h, zNear, zFar);
	//glUniformMatrix4fv(projection_matrix_loc, 1, GL_FALSE, &projection_matrix[0][0]);

	glm::mat4 view_matrix = 
		glm::lookAt(glm::vec3(camera_eye.X, camera_eye.Y, camera_eye.Z), 
					glm::vec3(camera_target.dX, camera_target.dY, camera_target.dZ),
					glm::vec3(camera_up.dX, camera_up.dY, camera_up.dZ));
	//glUniformMatrix4fv(view_matrix_loc, 1, GL_FALSE, &view_matrix[0][0]);

	glm::mat3 normal_matrix = glm::transpose(glm::inverse(glm::mat3(view_matrix)));
	//glUniformMatrix4fv(normal_matrix_loc, 1, GL_FALSE, &normal_matrix[0][0]);

	glUniformMatrix4fv(projection_matrix_loc, 1, GL_FALSE, &projection_matrix[0][0]);
	glUniformMatrix4fv(view_matrix_loc, 1, GL_FALSE, &view_matrix[0][0]);
	glUniformMatrix4fv(normal_matrix_loc, 1, GL_FALSE, &normal_matrix[0][0]);


	GLboolean save_color_mask[4];
	GLboolean save_depth_mask;
	glGetBooleanv(GL_COLOR_WRITEMASK, save_color_mask);
	glGetBooleanv(GL_DEPTH_WRITEMASK, &save_depth_mask);

	if (portalActivated[0]) {


		glm::mat4 portal_cam = portal_view(view_matrix, portals[0], portals[1]);
		glUniformMatrix4fv(view_matrix_loc, 1, GL_FALSE, &portal_cam[0][0]);
		glm::mat3 portal_normal_matrix = glm::transpose(glm::inverse(glm::mat3(portal_cam)));
		glUniformMatrix4fv(normal_matrix_loc, 1, GL_FALSE, &portal_normal_matrix[0][0]);
		level->displayScene(portal_cam);
		character->displayScene(portal_cam);
		//portals[0].displayScene(portal_cam);
		//portals[1].displayScene(portal_cam);

		//glFlush();


		glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
		glDepthMask(GL_TRUE);
		glStencilFunc(GL_NEVER, 0, 0xFF);
		glStencilOp(GL_INCR, GL_KEEP, GL_KEEP);  // draw 1s on test fail (always)
												 // draw stencil pattern
		glClear(GL_STENCIL_BUFFER_BIT);  // needs mask=0xFF
		glUniformMatrix4fv(view_matrix_loc, 1, GL_FALSE, &view_matrix[0][0]);
		portals[0].displayScene(view_matrix);

		glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
		glDepthMask(GL_TRUE);
		glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
		// Fill 1 or more 
		glStencilFunc(GL_LEQUAL, 1, 0xFF);

		glUniformMatrix4fv(view_matrix_loc, 1, GL_FALSE, &portal_cam[0][0]);
		// -Ready to draw main scene-

	}

	if (portalActivated[1]) {

		glm::mat4 portal_cam2 = portal_view(view_matrix, portals[1], portals[0]);
		glUniformMatrix4fv(view_matrix_loc, 1, GL_FALSE, &portal_cam2[0][0]);
		glm::mat3 portal_normal_matrix2 = glm::transpose(glm::inverse(glm::mat3(portal_cam2)));
		glUniformMatrix4fv(normal_matrix_loc, 1, GL_FALSE, &portal_normal_matrix2[0][0]);
		level->displayScene(portal_cam2);
		character->displayScene(portal_cam2);
		//portals[0].displayScene(portal_cam2);
		//portals[1].displayScene(portal_cam2);


		//glFlush();


		glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
		glDepthMask(GL_FALSE);
		glStencilFunc(GL_NEVER, 0, 0xFF);
		glStencilOp(GL_INCR, GL_KEEP, GL_KEEP);  // draw 1s on test fail (always)
												 // draw stencil pattern
		glClear(GL_STENCIL_BUFFER_BIT);  // needs mask=0xFF
		glUniformMatrix4fv(view_matrix_loc, 1, GL_FALSE, &view_matrix[0][0]);
		portals[1].displayScene(view_matrix);

		glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
		glDepthMask(GL_TRUE);
		glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
		// Fill 1 or more 
		glStencilFunc(GL_LEQUAL, 1, 0xFF);

		glUniformMatrix4fv(view_matrix_loc, 1, GL_FALSE, &portal_cam2[0][0]);
		// -Ready to draw main scene-

	}
	glColorMask(save_color_mask[0], save_color_mask[1], save_color_mask[2], save_color_mask[3]);
	glDepthMask(save_depth_mask);
/*
	if (portalActivated[0]) {


		//glFlush();

		GLboolean save_color_mask[4];
		GLboolean save_depth_mask;
		glGetBooleanv(GL_COLOR_WRITEMASK, save_color_mask);
		glGetBooleanv(GL_DEPTH_WRITEMASK, &save_depth_mask);

		glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
		glDepthMask(GL_FALSE);
		glStencilFunc(GL_NEVER, 0, 0xFF);
		glStencilOp(GL_INCR, GL_KEEP, GL_KEEP);  // draw 1s on test fail (always)
												 // draw stencil pattern
		glClear(GL_STENCIL_BUFFER_BIT);  // needs mask=0xFF
		glUniformMatrix4fv(view_matrix_loc, 1, GL_FALSE, &view_matrix[0][0]);
		portals[0].displayScene(view_matrix);

		glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
		glDepthMask(GL_TRUE);
		glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
		// Fill 1 or more //
		glStencilFunc(GL_LEQUAL, 1, 0xFF);

		glColorMask(save_color_mask[0], save_color_mask[1], save_color_mask[2], save_color_mask[3]);
		glDepthMask(save_depth_mask);

		glm::mat4 portal_cam = portal_view(view_matrix, portals[0], portals[1]);
		glUniformMatrix4fv(view_matrix_loc, 1, GL_FALSE, &portal_cam[0][0]);
		glm::mat3 portal_normal_matrix = glm::transpose(glm::inverse(glm::mat3(portal_cam)));
		glUniformMatrix4fv(normal_matrix_loc, 1, GL_FALSE, &portal_normal_matrix[0][0]);
		level->displayScene(portal_cam);
		character->displayScene(portal_cam);
		portals[0].displayScene(portal_cam);
		portals[1].displayScene(portal_cam);
		glUniformMatrix4fv(view_matrix_loc, 1, GL_FALSE, &portal_cam[0][0]);

		// -Ready to draw main scene-

	}


	if (portalActivated[1]) {

		//glFlush();

		GLboolean save_color_mask[4];
		GLboolean save_depth_mask;
		glGetBooleanv(GL_COLOR_WRITEMASK, save_color_mask);
		glGetBooleanv(GL_DEPTH_WRITEMASK, &save_depth_mask);

		glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
		glDepthMask(GL_FALSE);
		glStencilFunc(GL_NEVER, 0, 0xFF);
		glStencilOp(GL_INCR, GL_KEEP, GL_KEEP);  // draw 1s on test fail (always)
												 // draw stencil pattern
		glClear(GL_STENCIL_BUFFER_BIT);  // needs mask=0xFF
		glUniformMatrix4fv(view_matrix_loc, 1, GL_FALSE, &view_matrix[0][0]);
		portals[1].displayScene(view_matrix);

		glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
		glDepthMask(GL_TRUE);
		glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
		// Fill 1 or more //
		glStencilFunc(GL_LEQUAL, 1, 0xFF);

		glColorMask(save_color_mask[0], save_color_mask[1], save_color_mask[2], save_color_mask[3]);
		glDepthMask(save_depth_mask);

		glm::mat4 portal_cam2 = portal_view(view_matrix, portals[0], portals[1]);
		glUniformMatrix4fv(view_matrix_loc, 1, GL_FALSE, &portal_cam2[0][0]);
		glm::mat3 portal_normal_matrix2 = glm::transpose(glm::inverse(glm::mat3(portal_cam2)));
		glUniformMatrix4fv(normal_matrix_loc, 1, GL_FALSE, &portal_normal_matrix2[0][0]);
		level->displayScene(portal_cam2);
		character->displayScene(portal_cam2);
		portals[0].displayScene(portal_cam2);
		portals[1].displayScene(portal_cam2);
		glUniformMatrix4fv(view_matrix_loc, 1, GL_FALSE, &portal_cam2[0][0]);
		// -Ready to draw main scene-

	}*/

	glUniformMatrix4fv(projection_matrix_loc, 1, GL_FALSE, &projection_matrix[0][0]);
	glUniformMatrix4fv(view_matrix_loc, 1, GL_FALSE, &view_matrix[0][0]);
	glUniformMatrix4fv(normal_matrix_loc, 1, GL_FALSE, &normal_matrix[0][0]);

	// 3.6. Display the light
	// - point for point light
	// - segment for direction light
	glUniform1i(renderStyle_loc, 2);
	for (int i = 0; i<nbLight; i++) {
		lights[i].drawLight(shaderprogram1);
	}

	level->displayScene(view_matrix);
	character->displayScene(view_matrix);


	/*if (portalActivated[1]) {
		glClear(GL_DEPTH_BUFFER_BIT);
		GLboolean save_color_mask[4];
		GLboolean save_depth_mask;
		glGetBooleanv(GL_COLOR_WRITEMASK, save_color_mask);
		glGetBooleanv(GL_DEPTH_WRITEMASK, &save_depth_mask);
		glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
		glDepthMask(GL_TRUE);

		portals[1].displayScene(view_matrix);

		glColorMask(save_color_mask[0], save_color_mask[1], save_color_mask[2], save_color_mask[3]);
		glDepthMask(save_depth_mask);
	}*/

	//glFlush();
	glutSwapBuffers();
}

//This function is called from the main to initalize everything.
void init()
{
	
	camera_right = camera_forward.crossproduct(camera_up);
	camera_right.normalize();

	camera_target.dX = camera_eye.X + camera_forward.dX;
	camera_target.dY = camera_eye.Y + camera_forward.dY;
	camera_target.dZ = camera_eye.Z + camera_forward.dZ;

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
	
	// Portals
	/*double portal_vertices[] = {
		-1, -1, 0, 1,
		1, -1, 0, 1,
		-1,  1, 0, 1,
		1,  1, 0, 1,
	};
	for (unsigned int i = 0; i < sizeof(portal_vertices) / sizeof(portal_vertices[0]); i++) {
		portals[0].vertices.push_back(portal_vertices[i]);
		portals[1].vertices.push_back(portal_vertices[i]);
	}

	GLushort portal_elements[] = {
		0,1, 2, 3 , 4, 5, 6 ,7 ,8, 8, 9, 10, 11,4, 5 ,6 ,7 ,12, 13, 14, 15
	};
	for (unsigned int i = 0; i < sizeof(portal_elements) / sizeof(portal_elements[0]); i++) {
		portals[0].indices.push_back(portal_elements[i]);
		portals[1].indices.push_back(portal_elements[i]);
	}*/

	portals[0].myshaderprogram = shaderprogram1;
	portals[0].readScene("portal.obj");
	portals[0].createObjectBuffers();

	portals[1].myshaderprogram = shaderprogram1;
	portals[1].readScene("portal.obj");
	portals[1].createObjectBuffers();

	// Maya Scene
	level = new myObject3D(shaderprogram1);
	//level->readScene("test.obj");
	//level->readScene("TheCarnival.obj");
	level->readScene("portaltri4scale.obj");
	if (level->normals.size() == 0)
		level->computeNormals();
	level->createObjectBuffers();


	// Personnage iron man
	character = new myObject3D(shaderprogram1);
	character->readScene("IronMan2.obj");
	if (character->normals.size() == 0)
		character->computeNormals();
	character->createObjectBuffers();
	character->scale(0.008, 0.008, 0.008);
	character->translate(camera_eye.X, 2, camera_eye.Z);
	//character->rotate(0,1,0,180);

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
void animation(void)
{
	if (GLUTmouse[0] > Glut_w - 10) {
		teta -= 10 * sensitivity;
		//float rayon = sqrt(pow(dx, 2) + pow(dy, 2));
		r_temp = cos(phi * PI / 180);
		camera_forward.dX = sin(teta * PI / 180.0) * r_temp;
		camera_forward.dY = sin(phi * PI / 180.0);
		camera_forward.dZ = cos(teta * PI / 180.0) * r_temp;

		if (character != nullptr) {
			character->translate(-camera_eye.X, 0, -camera_eye.Z);
			character->rotate(0, 1, 0, -arm_angle);
			character->rotate(0, 1, 0, teta);
			character->translate(camera_eye.X, 0, camera_eye.Z);
			arm_angle = teta;
		}

		//double theta = 4.0 * (fabs(vx) + fabs(vy));

		camera_right = camera_forward.crossproduct(camera_up);
		camera_right.normalize();

		camera_target.dX = camera_eye.X + camera_forward.dX;
		camera_target.dY = camera_eye.Y + camera_forward.dY;
		camera_target.dZ = camera_eye.Z + camera_forward.dZ;
	}
	else if (GLUTmouse[0] < 15) {
		teta += 10 * sensitivity;
		//float rayon = sqrt(pow(dx, 2) + pow(dy, 2));
		r_temp = cos(phi * PI / 180);
		camera_forward.dX = sin(teta * PI / 180.0) * r_temp;
		camera_forward.dY = sin(phi * PI / 180.0);
		camera_forward.dZ = cos(teta * PI / 180.0) * r_temp;

		if (character != nullptr) {
			character->translate(-camera_eye.X, 0, -camera_eye.Z);
			character->rotate(0, 1, 0, -arm_angle);
			character->rotate(0, 1, 0, teta);
			character->translate(camera_eye.X, 0, camera_eye.Z);
			arm_angle = teta;
		}

		//double theta = 4.0 * (fabs(vx) + fabs(vy));

		camera_right = camera_forward.crossproduct(camera_up);
		camera_right.normalize();

		camera_target.dX = camera_eye.X + camera_forward.dX;
		camera_target.dY = camera_eye.Y + camera_forward.dY;
		camera_target.dZ = camera_eye.Z + camera_forward.dZ;
	}



	if (portalsDirectionIteraction[0] != 0) {
		//portals[0].translate(portalsDirection[0].dX, portalsDirection[0].dY, portalsDirection[0].dZ);
		portals[0].translate(portalsDirection[0].dX, 0, portalsDirection[0].dZ);
		portalsDirectionIteraction[0]--;
	}
	if (portalsDirectionIteraction[1] != 0) {
		//portals[1].translate(portalsDirection[1].dX, portalsDirection[1].dY, portalsDirection[1].dZ);
		portals[1].translate(portalsDirection[1].dX, 0, portalsDirection[1].dZ);
		portalsDirectionIteraction[1]--;
	}

	glutPostRedisplay();
}


int main(int argc, char* argv[]) {
	glutInit(&argc, argv);
	//glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA | GLUT_DEPTH );
	glutInitContextVersion(2, 0);
	glutInitDisplayMode(GLUT_RGBA | GLUT_ALPHA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_STENCIL);
	glutCreateWindow("My OpenGL Application");
	   
	glewInit() ; 
	glutReshapeWindow(Glut_w, Glut_h);
	
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(keyboard2);
	glutPassiveMotionFunc(mousedrag) ;
	glutMouseFunc(mouse) ;
	glutMouseWheelFunc(mouseWheel);
	glutIdleFunc(animation);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_STENCIL_TEST);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_CULL_FACE);
	glDepthFunc(GL_LESS) ;

	//glEnable(GL_POLYGON_OFFSET_FILL);
	//glPolygonOffset(1, 0);
	
	init();

	glutMainLoop();
	return 0;
}