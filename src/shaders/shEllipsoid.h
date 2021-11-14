#ifndef __Ellipsoid__sh__
#define __Ellipsoid__sh__

inline std::string Ellipsoid_fsh() 
{
	std::string str = 
	"#version 330 core\n"\
	"in vec4 vColor;\n"\
	"in vec2 vTex;\n"\
	"in vec4 vPos;\n"\
	"in vec4 vExtra;\n"\
	"\n"\
	"uniform sampler2D pic_tex;\n"\
	"\n"\
	"out vec4 fragColor;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"	fragColor = vColor;\n"\
	"	vec2 screenPos = (vPos.xy / vPos.w + 1) * 100;\n"\
	"	ivec2 rounded = ivec2(screenPos);\n"\
	"	if (vExtra.x > 0.5 && (rounded.x + rounded.y) % 2 == 0)\n"\
	"	{\n"\
	"		fragColor.xyz = vec3(1., 1., 1.);\n"\
	"	}\n"\
	"\n"\
	"\n"\
	"}\n";
	return str;
}

inline std::string Ellipsoid_vsh()
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
	"out vec4 vExtra;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"    vec4 pos = vec4(position[0], position[1], position[2], 1.0);\n"\
	"    vPos = projection * model * pos;\n"\
	"    gl_Position = vPos;\n"\
	"	 vColor = color;\n"\
	"	 vExtra = extra;\n"\
	"}";
	return str;
}

#endif
