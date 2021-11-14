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

#ifndef __Tent__Tent__
#define __Tent__Tent__

#include <h3dsrc/SlipObject.h>

class Workbench;

class Tent : public SlipObject
{
public:
	Tent(Workbench *bench);
	
	void makeBase(std::string strain);

	void setWorkbench(Workbench *bench)
	{
		_bench = bench;
	}
	
	Workbench *bench()
	{
		return _bench;
	}

private:
	typedef struct
	{
		vec3 p;
		double l;
	} Distance;

	static bool compare_distance(Distance &a, Distance &b)
	{
		return (a.l < b.l);
	}

	std::vector<vec3> _ps;
	void closestPoints(vec3 v, std::vector<Distance> *results);

	Workbench *_bench;

};

#endif
