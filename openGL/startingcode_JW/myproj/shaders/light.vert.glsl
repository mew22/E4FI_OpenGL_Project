#version 330 core

layout(location = 0) in vec4 vertex_modelspace;
layout(location = 2) in vec3 normal_modelspace;

out vec3 mynormal;
out vec4 myvertex;

uniform mat4 myprojection_matrix;
uniform mat4 myview_matrix;
uniform mat3 mynormal_matrix;
uniform mat4 mymodel_matrix;

void main() {

	 //gl_Position = myprojection_matrix * myview_matrix * vertex_modelspace; 
	 gl_Position = myprojection_matrix * myview_matrix * mymodel_matrix * vertex_modelspace;

	mynormal = normal_modelspace;
	myvertex = vertex_modelspace;
}
