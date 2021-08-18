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
#include <hcsrc/maths.h>
#include <h3dsrc/Icosahedron.h>
#include <c4xsrc/GLAxis.h>
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

typedef struct
{
	Mutation *mut;
	double val;
} OrderMut;

bool order_by_val(const OrderMut &a, const OrderMut &b)
{
	return (a.val < b.val);
}

void Scatter::reorderMutations(std::vector<Mutation *> &muts)
{
	std::vector<OrderMut> sorting;
	OrderMut o;

	for (size_t i = 0; i < muts.size(); i++)
	{
		muts[i]->calculateCloud();
		vec3 c = muts[i]->centre();
		o.mut = muts[i];
		o.val = vec3_length(c);
		sorting.push_back(o);
	}

	std::sort(sorting.begin(), sorting.end(), order_by_val);
	
	for (size_t i = 0; i < muts.size(); i++)
	{
		muts[i] = sorting[i].mut;
	}
}

Icosahedron *Scatter::ico(mat3x3 tensor, vec3 pos, double colour, int tri)
{
	double r, g, b;
	val_to_cluster4x_colour(colour, &r, &g, &b);
	r /= 350; g /= 350; b /= 350;

	Icosahedron *m = new Icosahedron();
	for (size_t i = 0; i < tri; i++)
	{
		m->triangulate();
	}
	m->setColour(0.5, 0.5, 0.5);
	m->setAlpha(0.3);
	m->rotateByMatrix(tensor);
	m->addToVertices(pos);
	double mute_r = 0.6 + (r - 0.5) * 0.6;
	double mute_g = 0.6 + (g - 0.5) * 0.6;
	double mute_b = 0.6 + (b - 0.5) * 0.6;
	m->setColour(mute_r, mute_g, mute_b);
	m->setAlpha(0.3);

	return m;
}

void Scatter::strainWalk(Strain *strain)
{
	vec3 pos = empty_vec3();
	GLAxis *walk = new GLAxis(make_vec3(1, 0, 0));
	walk->clearVertices();
	_objs.push_back(walk);
	walk->addVertex(empty_vec3());
	
	std::vector<Mutation *> list = strain->list();
	reorderMutations(list);

	for (size_t n = 0; n < list.size(); n++)
	{
		if (list[n]->silenced())
		{
			continue;
		}
		
		vec3 c = list[n]->centre();
		pos += c;

		if (list[n]->trials() >= 3)
		{
			mat3x3 t = list[n]->tensor();
			double c = list[n]->averageEase() * 4;
			Icosahedron *m = ico(t, pos, c);
			_objs.push_back(m);
		}
		
		std::string str = list[n]->str();
		
		if (n == list.size() - 1)
		{
			str += " (" + strain->name() + ")";
		}

		addPoint(pos, str);

		walk->addIndex(-2);
		walk->addIndex(-1);

		walk->addVertex(pos);
	}

	if (walk->vertexCount() > 2)
	{
		walk->addIndex(-2);
		walk->addIndex(-1);
	}
}

void Scatter::clearBalls()
{
	// get rid of old meshes
	lockMutex();
	for (size_t i = 0; i < _objs.size(); i++)
	{
		delete _objs[i];
	}
	_objs.clear();
	unlockMutex();

}

void Scatter::populateFromMutations(const std::vector<Mutation *> &muts)
{
	repopulate();
	clearBalls();
	_renderType = GL_POINTS;

	for (size_t n = 0; n < muts.size(); n++)
	{
		if (muts[n]->silenced())
		{
			continue;
		}
		
		muts[n]->calculateCloud();
		vec3 c = muts[n]->centre();

		if (muts[n]->trials() >= 3)
		{
			mat3x3 t = muts[n]->tensor();
			double col = muts[n]->averageEase() * 4;
			Icosahedron *m = ico(t, c, col, 2);
			_objs.push_back(m);
		}

		addPoint(c, muts[n]->str());

		continue;

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
	clearBalls();
	
	if (strains.size() == 1)
	{
		strainWalk(strains[0]);
		return;
	}

	_renderType = GL_POINTS;
	for (size_t n = 0; n < strains.size(); n++)
	{
		vec3 point = strains[n]->aveDirection();
		double ease = strains[n]->averageEase() * 4;
		
		if (strains[n]->trials() >= 3)
		{
			strains[n]->calculateCloud();
			mat3x3 tensor = strains[n]->tensor();
			Icosahedron *m = ico(tensor, point, ease, 2);
			_objs.push_back(m);
		}

		addPoint(point, strains[n]->name());
		int idx = vertexCount() - 1;
		Helen3D::Vertex *v = &_vertices[idx];
		/*
		v->color[0] = r;
		v->color[1] = g;
		v->color[2] = b;

		_texts[idx]->setColour(r, g, b, 1);
		*/
		_texts[idx]->prepare();
	}
}

void Scatter::render(SlipGL *gl)
{
	Plot3D::render(gl);

	lockMutex();
	for (size_t i = 0; i < _objs.size(); i++)
	{
		_objs[i]->render(gl);
		_objs[i]->reorderIndices();
	}
	unlockMutex();
}
