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

#include "Settings.h"
#include "Workbench.h"
#include "Strain.h"
#include <hcsrc/FileReader.h>

Settings::Settings(Workbench *workbench, std::string filename)
{
	_filename = filename;
	_valid = file_exists(_filename);
	_bench = workbench;
	
	if (_valid)
	{
		try
		{
			_contents = get_file_contents(_filename);
		}
		catch (int e)
		{
			_valid = false;
		}
	}
	
	std::cout << "Valid: " << _valid << std::endl;
}

void Settings::apply()
{
	std::vector<std::string> lines = split(_contents, '\n');
	std::cout << "Applying" << std::endl;

	for (size_t i = 0; i < lines.size(); i++)
	{
		std::string line = lines[i];
		std::vector<std::string> bits = split(line, ' ');
		
		if (bits.size() <= 0)
		{
			continue;
		}
		
		if (bits[0] == "set")
		{
			if (bits.size() < 3)
			{
				continue;
			}
			std::string keyword = bits[1];
			std::string value = bits[2];
			_dict[keyword] = value;
		}
		
		Strain *str = _bench->strain(bits[0]);
		if (str == NULL)
		{
			std::cout << "Could not find " << bits[0] << std::endl;
			continue;
		}
		
		if (bits.size() < 4)
		{
			continue;
		}

		double red = atof(bits[1].c_str());
		double green = atof(bits[2].c_str());
		double blue = atof(bits[3].c_str());
		
		str->setColour(red, green, blue);
		std::cout << "Setting " << str->name() << " colour" << std::endl;

		str->setChequered(false);
		
		if (bits.size() < 5)
		{
			continue;
		}
		
		int next = 4;
		
		if (bits[next] == "#")
		{
			std::cout << "Setting " << str->name() << " to chequered" << std::endl;
			str->setChequered(true);
			next++;
		}
		
		if (bits.size() <= next)
		{
			continue;
		}
		
		std::string display_name = bits[next];
		str->setDisplayName(display_name);
		next++;
	}

}
