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

#include "Scatter.h"
#include "Mutation.h"
#include "Strain.h"
#include <c4xsrc/shaders/Blob_fsh.h>
#include <c4xsrc/shaders/Blob_vsh.h>

Scatter::Scatter() : Plot3D()
{
	setName("Scatter");
	_renderType = GL_POINTS;
	_fString = pointFsh();
	_vString = pointVsh();
	setShowText(true);
	setPointSize(5);
}

void Scatter::populateFromMutations(const std::vector<Mutation *> &muts)
{
	repopulate();

	for (size_t n = 0; n < muts.size(); n++)
	{
		if (muts[n]->silenced())
		{
			continue;
		}

		for (size_t i = 0; i < muts[n]->trials(); i++)
		{
			std::vector<double> results = muts[n]->savedVector(i);
			
			if (results.size() < 3)
			{
				results.resize(3);
			}
			
			vec3 point = make_vec3(results[_a], results[_b], results[_c]);
			addPoint(point, muts[n]->str());
		}
	}
}

void Scatter::populateFromStrains(const std::vector<Strain *> &strains)
{
	repopulate();

	for (size_t n = 0; n < strains.size(); n++)
	{
		std::vector<double> results = strains[n]->direction();

		if (results.size() < 3)
		{
			results.resize(3);
		}

		vec3 point = make_vec3(results[_a], results[_b], results[_c]);
		addPoint(point, strains[n]->name());
	}
}
