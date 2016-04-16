#version 330 core

uniform int myrenderStyle;

in vec3 mynormal;
in vec4 myvertex;

uniform mat4 myview_matrix;
uniform mat3 mynormal_matrix;
uniform mat4 mymodel_matrix;

uniform vec4 mylight_position = vec4(0,2,2,1);
uniform vec4 mylight_color = vec4(1,1,1,0); 

uniform vec4 mymodel_kd;
uniform vec4 mymodel_ks;
uniform vec4 mymodel_ka;
uniform float mymodel_ns;

uniform vec4 light_colors[16];
uniform vec4 light_position[16];
uniform vec4 light_direction[16];
uniform int  light_type[16];



void main (void)
{   
	float ang =	1;
	float spotAngle = 22.5;
	vec3 eyepos = vec3(0,0,0);
	int puissance = 30;
	vec4 _mypos = myview_matrix * myvertex;
	vec3 mypos = _mypos.xyz / _mypos.w;

	vec3 normal = normalize( mynormal_matrix * mynormal );

	vec3 mypos_to_lightpos;

	for(int i = 0;i<16;i++){
		vec3 direction = normalize(mynormal_matrix*light_direction[i].xyz);
		vec4 _lightpos = myview_matrix*light_position[i];
		vec3 lightpos = _lightpos.xyz/_lightpos.w;
		if( light_type[i]== 1){
			mypos_to_lightpos = normalize( lightpos - mypos );
		}else if( light_type[i]== 2){
			mypos_to_lightpos = direction;
		}else if( light_type[i]== 3) {
			mypos_to_lightpos = normalize( lightpos - mypos );

			vec3 posTemp = normalize( mypos - lightpos );
			ang = pow(dot(direction,posTemp),puissance);
		}


		vec3 mypos_to_eyepos = normalize( eyepos - mypos );

	
		gl_FragColor += light_colors[i] * ang * mymodel_kd * dot( normal, mypos_to_lightpos );

		vec3 reflectedray = reflect( -mypos_to_lightpos, normal );

		//specular color
		gl_FragColor += light_colors[i] * mymodel_ks * ang *
					   pow( max( dot( reflectedray, mypos_to_eyepos ), 0.0), 200) ;

	
	}

}
