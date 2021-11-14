#ifndef __Serum_contour__
#define __Serum_contour__

inline std::string Contour_vsh()
{
	std::string str = 
	"#version 330 core\n"\
	"in vec3 normal;\n"\
	"in vec3 position;\n"\
	"in vec4 color;\n"\
	"in vec4 extra;\n"\
	"in vec2 tex;\n"\
	"\n"\
	"uniform mat4 model;\n"\
	"uniform mat4 projection;\n"\
	"uniform float time;\n"\
	"\n"\
	"out vec4 vColor;\n"\
	"out vec4 vPos;\n"\
	"out vec2 vTex;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"    vec4 pos = vec4(position[0], position[1], position[2], 1.0);\n"\
	"	 vec4 norm4 = vec4(0., 0., 1.0, 0.0);\n"\
	"	 pos = model * pos;\n"\
	"	 vec4 lightpos = vec4(pos[0], pos[1], pos[2], 1);\n"\
	"    float mag = abs(dot(normalize(norm4), normalize(lightpos)));"\
	"    vPos = vec4(position, 1.0);\n"\
	"    gl_Position = projection * pos;\n"\
	"	 vColor = vec4(color.xyz, mag);\n"\
	"    vTex = tex;\n"\
	"}";
	return str;
}

inline std::string Contour_fsh() 
{
	std::string str = 
	"#version 330 core\n"\
	"in vec4 vColor;\n"\
	"in vec2 vTex;\n"\
	"in vec4 vPos;\n"\
	"\n"\
	"uniform sampler2D pic_tex;\n"\
	"\n"\
	"out vec4 fragColor;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"	float z = vPos.z / 2;\n"\
	"	if (z >= 2)\n"\
	"	{\n"\
	"		fragColor = vec4(1., 1., z - 2, 1.);\n"\
	"	}\n"\
	"	else if (z >= 1)\n"\
	"	{\n"\
	"		fragColor = vec4(1., z - 1., 0.0, 1);\n"\
	"	}\n"\
	"	else if (z >= 0)\n"\
	"	{\n"\
	"		fragColor = vec4(z, 0.0, 0.0, 1);\n"\
	"	}\n"\
	"\n"\
	"\n"\
	"\n"\
	"}\n";
	return str;
}


#endif
