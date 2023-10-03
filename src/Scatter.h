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

#ifndef __serum__Scatter__
#define __serum__Scatter__

#include <h3dsrc/Plot3D.h>
#include <hcsrc/FileReader.h>

class Mutation;
class Strain;
class Icosahedron;

class Scatter : public Plot3D
{
public:
	Scatter();
	
	void populate() {};

	void populateFromMutations(const std::vector<Mutation *> &muts,
	                           float scale);
	void populateFromStrains(const std::vector<Strain *> &strains,
	                         float scale);
	void strainWalk(Strain *strain);
	
	virtual size_t axisCount()
	{
		return 3;
	}
	
	virtual std::string axisLabel(int i)
	{
		return "axis_" + i_to_str(i);
	}
	
	virtual void render(SlipGL *gl);
private:
	Icosahedron *ico(mat3x3 tensor, vec3 pos, double colour = 0, int tri = 1);
	void reorderMutations(std::vector<Mutation *> &muts);
	void clearBalls();
	std::vector<SlipObject *> _objs;

};

#endif
