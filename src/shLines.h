#ifndef __Serum__lines__
#define __Serum__lines__

inline std::string Fade_vsh()
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
	"    vec4 pos = vec4(position.xyz, 1.0);\n"\
	"	 mat3 rot = mat3(model);\n"\
	"	 vec4 norm4 = vec4(rot * normal.xyz, 0.0);\n"\
	"	 vec4 dir = vec4(0., 0., 1., 0.);\n"\
	"	 pos = model * pos;\n"\
	"    float mag = (dot(normalize(norm4), dir));"\
	"    vPos = pos;\n"\
	"    gl_Position = projection * pos;\n"\
	"	 vColor = vec4(color.xyz * mag, 1.);\n"\
	"    vTex = tex;\n"\
	"}";
	return str;
}

inline std::string Fade_fsh() 
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
	"	float z = vPos.z;\n"\
	"	z = -(z + 60) / 10;\n"\
	"	z = max(0, z);\n"\
	"	z = min(1, z);\n"\
	"   vec3 toWhite = vec3(1., 1., 1.) - vColor.xyz;\n"\
	"	fragColor = vec4(vColor.xyz + toWhite * (z), 1.);\n"\
	"\n"\
	"\n"\
	"\n"\
	"}\n";
	return str;
}

#endif
