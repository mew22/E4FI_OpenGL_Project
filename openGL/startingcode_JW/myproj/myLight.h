class myLight
{
	public:
	GLfloat position[4];
	GLfloat color[4];
	GLfloat direction[3];
	GLfloat type; // 0 (point), 1 (directional), 2 (spot)


	void drawLight(GLuint shader)
	{
		float length = sqrt(direction[0] * direction[0] + direction[1] * direction[1]+ direction[2] * direction[2]);
		if (length != 0) for (int i=0;i<3;i++) direction[i] /= length;
		glPointSize(6.0);
		glUniform4fv(glGetUniformLocation(shader, "kd"), 1, &color[0]);
		if (type == 0 || type == 2) {
			glBegin(GL_POINTS);
			glVertex3f(position[0], position[1], position[2]);
			glEnd();
		}
		if (type == 1){
			glBegin(GL_LINES);
			glVertex3f(0,0,0);
			glVertex3f(direction[0] , direction[1] , direction[2] );
			glEnd();
		}
		if (type == 2){
			glBegin(GL_LINES);
			glVertex3f(position[0], position[1], position[2]);
			glVertex3f(position[0] + direction[0] / 2.0, position[1] + direction[1] / 2.0,
			position[2] + direction[2] / 2.0);
			glEnd();
		}
	}


};


