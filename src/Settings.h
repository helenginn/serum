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

#ifndef __serum__Settings__
#define __serum__Settings__

#include <string>
#include <map>

class Workbench;

class Settings
{
public:
	Settings(Workbench *workbench, std::string filename);
	
	bool hasValue(const std::string &key)
	{
		return _dict.count(key);
	}
	
	const std::string &valueFor(const std::string &key) const
	{
		return _dict.at(key);
	}

	bool isValid()
	{
		return _valid;
	}
	
	void apply();
private:
	Workbench *_bench;
	std::string _filename;
	std::string _contents;
	bool _valid;

	std::map<std::string, std::string> _dict; 
};

#endif
