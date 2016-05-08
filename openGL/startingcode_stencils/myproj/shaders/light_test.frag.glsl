#version 330 core

uniform int myrenderStyle;

in vec3 mynormal;
in vec4 myvertex;
in vec2 mytextcoord;

uniform mat4 myview_matrix;
uniform mat3 mynormal_matrix;
uniform mat4 mymodel_matrix;

uniform vec4 mymodel_kd;
uniform vec4 mymodel_ks;
uniform vec4 mymodel_ka;
uniform float mymodel_ns;

// Création de plusieurs lights
uniform vec4 light_colors[16];
uniform vec4 light_position[16];
uniform vec3 light_direction[16];
uniform int light_type[16];
uniform int nbLights;

// Avant, pour afficher une lumière
vec4 mylight_position = vec4(0,2,2,1);
vec4 mylight_color = vec4(1,1,1,0); 
vec4 kd = vec4(1,1,1,0);
vec4 ks = vec4(1,1,1,0);

uniform sampler2D tex;
uniform int using_textures;


void main (void)
{   
	// replace Kd by texture()
	// pass boolean to know if using texture
	if(using_textures == 1)
		mymodel_kd = texture(tex, mytextcoord.st);
	vec3 eyepos = vec3(0,0,0);

	vec4 _mypos = myview_matrix * mymodel_matrix * myvertex;
	vec3 mypos = _mypos.xyz / _mypos.w;

	vec3 normal = normalize( mynormal_matrix * mynormal );


	if (myrenderStyle == 0) gl_FragColor = vec4(1,0,0,0);

	// Silhouete 
	if (myrenderStyle == 1) {
	
		vec4 final_pos = myview_matrix * myvertex;

		vec3 really_final_vertex = final_pos.xyz/final_pos.w;
		vec3 final_normal = mynormal_matrix*normal;

		if(dot(really_final_vertex,final_normal) < 0.1) 
			gl_FragColor = vec4(0,0,1,0);
	}

	if (myrenderStyle == 2)		gl_FragColor = mymodel_kd;
	
	
	// 3.5
	vec3 mypos_to_lightpos;
	vec3 mypos_to_eyepos = normalize( eyepos - mypos );

	gl_FragColor = vec4(0,0,0,0);
	for (int i = 0; i<nbLights; i++){
		vec3 lightpos = light_position[i].xyz ;
		vec3 lightdir = light_direction[i].xyz;
		
		// Point Light
		if(light_type[i] == 0){
			mypos_to_lightpos = normalize( lightpos - mypos );


			//diffuse color
			gl_FragColor += light_colors[i] * mymodel_kd * dot( normal, mypos_to_lightpos );

		
			vec3 reflectedray = reflect( -mypos_to_lightpos, normal );
			//specular color
			gl_FragColor += light_colors[i] * mymodel_ks * 
					pow( max( dot( reflectedray, mypos_to_eyepos ), 0.0), 200) ;
		}

		
		// Directionnal light
		else if(light_type[i] == 1){
			mypos_to_lightpos = normalize( mynormal_matrix* lightdir  ); // The directionnal light must be multiplied with the normal matrix to be appied on only 1 direction of the apple



			//diffuse color
			gl_FragColor += light_colors[i] * mymodel_kd * dot( normal, mypos_to_lightpos );

		
			vec3 reflectedray = reflect( -mypos_to_lightpos, normal );
			//specular color
			gl_FragColor += light_colors[i] * mymodel_ks * 
					pow( max( dot( reflectedray, mypos_to_eyepos ), 0.0), 200) ;

		}

		// Spotlight
		else if(light_type[i] == 2){
			// We have mypos_to_lightpos = p*cos(theta)
			float cos_theta = dot(lightdir,mypos-lightpos); // mypos-lightpos : vector from lightpos to mypos

			mypos_to_lightpos = normalize( lightpos*cos_theta );

			
			//diffuse color
			gl_FragColor += light_colors[i] * mymodel_kd * dot( normal, mypos_to_lightpos );

		
			vec3 reflectedray = reflect( -mypos_to_lightpos, normal );
			//specular color
			gl_FragColor += light_colors[i] * mymodel_ks * 
					pow( max( dot( reflectedray, mypos_to_eyepos ), 0.0), 200) ;			
		}
		
		
	}

 	

}