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

inline std::string Foggy_fsh() 
{
	std::string str = 
	"#version 330 core\n"\
	"in vec4 vColor;\n"\
	"in vec4 vPos;\n"\
	"in vec3 mPos;\n"\
	"in vec4 vExtra;\n"\
	"\n"\
	"uniform sampler2D pic_tex;\n"\
	"\n"\
	"out vec4 fragColor;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"	fragColor = vColor;\n"\
	"	\n"\
	"	\n"\
	"	float amount = min(1 - ((mPos.z + 1) / 2), 1);\n"\
	"	fragColor.xyz += (1 - fragColor.xyz) * amount;\n"\
	"\n"\
	"\n"\
	"}\n";
	return str;
}

inline std::string Foggy_vsh()
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
	"out vec3 mPos;\n"\
	"out vec4 vExtra;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"    vec4 pos = vec4(position[0], position[1], position[2], 1.0);\n"\
	"	 vec4 norm4 = vec4(normal[0], normal[1], normal[2], 1.0);\n"\
	"	 pos = model * pos;\n"\
	"	 vec3 norm3 = mat3(model) * norm4.xyz;\n"\
	"	 vec3 lightpos = vec3(pos[0], pos[1], pos[2]);\n"\
	"    float mag = dot(normalize(norm3), normalize(lightpos));"\
	"    vPos = projection * pos;\n"\
	"    mPos = mat3(model) * position;\n"\
	"    gl_Position = vPos;\n"\
	"	 vColor = vec4(color.xyz * (1 - mag), color[3]);\n"\
	"	 vExtra = extra;\n"\
	"}";
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
