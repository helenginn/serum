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

#ifndef __serum__Serum__
#define __serum__Serum__

//#include <vector>
#include <string>

class Strain;

class Serum
{
public:
	Serum(std::string name, Strain *strain);

	const std::string &name()
	{
		return _name;
	}
	
	Strain *strain() 
	{
		return _strain;
	}
	
	double strength()
	{
		return _strength;
	}
	
	double *strengthPtr()
	{
		return &_strength;
	}
	
	double offset()
	{
		return _offset;
	}
	
	double *offsetPtr()
	{
		return &_offset;
	}
private:
	std::string _name;
	Strain *_strain;
	double _strength;
	double _offset;

};

#endif
