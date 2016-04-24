#version 330 core

uniform int myrenderStyle;

in vec3 mynormal;
in vec4 myvertex;
in vec2 mytextcoord;
in vec4 mycolor;

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

uniform sampler2D tex;
uniform int using_textures;


vec4 computeLight(  const in vec3 eyepos, const in vec3 mypos, const in vec3 lightpos, const in vec3 mylight_direction, const in int mylight_type,
					in vec3 normal, const in vec4 ka, const in vec4 kd, 
					const in vec4 ks, const in float s )
{
	vec3 eyedir = normalize(eyepos - mypos) ;
	vec3 lightdir = normalize (lightpos - mypos) ;

	if (dot(normal, eyedir)<0) normal = -normal;
	vec3 reflectdir = normalize( reflect(-lightdir, normal) );
	vec4 color = mycolor;
	
    //point-lights
	if (mylight_type == 1) {
	      color =  ka + kd * mylight_color * max( dot(lightdir, normal), 0.0) +
				   ks * mylight_color * pow( max( dot(reflectdir, eyedir), 0.0 ), s );
    } 
	//directinal lights
	else if (mylight_type == 2)
	{
		color =  ka + kd * mylight_color * max( dot(-mynormal_matrix*mylight_direction, normal), 0.0) +
		         ks * mylight_color * pow( max( dot(normalize( reflect(mynormal_matrix*mylight_direction, normal) ), eyedir), 0.0 ), s );
		//the above is correct, as have to mult by normalmatrix.however, if want fixed direction:
	    //gl_FragColor =  gl_FrontMaterial.ambient + 
		//			gl_FrontMaterial.diffuse * mylight_color * max( dot(-mylight_direction, normal), 0.0) +
		//			gl_FrontMaterial.specular * mylight_color * pow( max( dot(normalize( reflect(mylight_direction, normal) ), eyedir), 0.0 ), gl_FrontMaterial.shininess );
	} 
	//spot-lights
	else if (mylight_type == 3)
	{
	    color =  ka + kd * mylight_color * max( dot(lightdir, normal), 0.0) +
			     ks * mylight_color * pow( max( dot(reflectdir, eyedir), 0.0 ), s );
		color = gl_FragColor*pow(max(dot(-lightdir,normalize(mylight_direction)),0.0), 2000);
	}

	return color;
}





void main (void)
{   
    vec4 kd = mymodel_kd.xyzw;

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

	if (myrenderStyle == 2)		gl_FragColor = kd;
	
	
	// 3.5
	vec3 mypos_to_lightpos;
	vec3 mypos_to_eyepos = normalize( eyepos - mypos );

	// Textured
	if(using_textures == 1){
		kd = texture(tex, mytextcoord.st);
	}

	gl_FragColor = vec4(0,0,0,0);
	for (int i = 0; i<nbLights; i++){
		vec3 lightpos = light_position[i].xyz ;
		vec3 lightdir = light_direction[i].xyz;
		gl_FragColor += computeLight(eyepos, mypos, lightpos, lightdir, light_type[i],normal, mymodel_ka, kd, mymodel_ks, mymodel_ns);
	}

}

