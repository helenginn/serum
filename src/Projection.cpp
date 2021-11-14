// serum
// Copyright (C) 2019 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.

#include "Projection.h"
#include <hcsrc/FileReader.h>
#include <h3dsrc/SlipGL.h>

Projection::Projection(std::string filename)
{
	std::string details;
	try
	{
		details = get_file_contents(filename);
	}
	catch (int e)
	{
		std::cout << "Could not open file, " << filename << std::endl;
		return;
	}
	
	std::vector<std::string> components = split(details, ' ');

	mat4x4 m = make_mat4x4();

	for (size_t i = 0; i < components.size() && i < 16; i++)
	{
		double d = atof(components[i].c_str());
		m.vals[i] = d;
	}
	
	_mat = m;
	_centre = empty_vec3();
	
	if (components.size() >= 19)
	{
		_centre.x = atof(components[16].c_str());
		_centre.y = atof(components[17].c_str());
		_centre.z = atof(components[18].c_str());
	}
}

void Projection::applyToGL(SlipGL *gl)
{
	gl->setModel(_mat);
	gl->changeCentre(_centre);
}
