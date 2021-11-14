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

#include "Tent.h"
#include "Workbench.h"
#include "Strain.h"
#include "Serum.h"
#include <h3dsrc/Icosahedron.h>
#include <h3dsrc/Mesh.h>
#include <h3dsrc/shaders/vStructure.h>
#include <h3dsrc/shaders/fStructure.h>
#include "shContours.h"
#include <iostream>

Tent::Tent(Workbench *bench) : SlipObject()
{
	_vString = Structure_vsh();
	_fString = Structure_fsh();
	_bench = bench;
}

void Tent::makeBase(std::string str)
{
	Icosahedron *ico = new Icosahedron();
	ico->triangulate();
	ico->triangulate();
	ico->triangulate();
	ico->setColour(0.0, 0.0, 0.0);
	mat3x3 scale = make_mat3x3();
	mat3x3_scale(&scale, 4, 4, 0.2);
	ico->rotateByMatrix(scale);

//	appendObject(ico);
	delete ico;

	Strain *strain = _bench->strain(str);
	std::cout << strain->name() << std::endl;

	for (size_t i = 0; i < _bench->strainCount(); i++)
	{
		Strain *challenge = _bench->strain(i);
		
		if (challenge->childrenCount() == 0)
		{
			continue;
		}

		double ave = 0;
		double count = 0;
		for (size_t k = 0; k < strain->childrenCount(); k++)
		{
			Serum *test = strain->child(k);
			double val = challenge->serumValue(test);
			
			if (val != val)
			{
				continue;
			}

			ave += val;
			count++;
		}
		
		ave /= count;

		if (ave <= 0.01 || ave != ave)
		{
			continue;
		}

		std::cout << std::endl;

		vec3 point = challenge->aveDirection();
		vec3_mult(&point, 3);
		mat3x3 scale = make_mat3x3();
		double height = ave;
		point.z = 0;
		mat3x3_scale(&scale, 0.2, 0.2, 1);
		Icosahedron *ico = new Icosahedron();
		ico->triangulate();

		ico->rotateByMatrix(scale);
		ico->setPosition(point);
		
		for (size_t j = 0; j < ico->vertexCount(); j++)
		{
			bool up = (ico->vertex(j).pos[2] > 0);
			ico->vPointer()[j].pos[2] = (up ? height * (2.2) : 0);
			ico->vPointer()[j].pos[2] = (up ? 1.0 : 0);
		}
		point.z = height * 2;
		_ps.push_back(point);

		appendObject(ico);
		delete ico;
	}
	
	double r = 3.8;
	for (double a = 0; a <= 361; a += 10)
	{
		double rad = M_PI * a / 180.;
		double x = r * cos(rad);
		double y = r * sin(rad);
		vec3 v = make_vec3(x, y, 0);
		_ps.push_back(v);
	}

	Mesh *m = makeMesh(6);

	{
		mat3x3 scale = make_mat3x3();
		mat3x3_scale(&scale, 2, 2, 1);
		m->rotateByMatrix(scale);
	}
	
	for (size_t i = 0; i < m->indexCount(); i+=3)
	{
		bool kill = false;
		for (size_t j = 0; j < 3; j++)
		{
			Helen3D::Vertex vert = m->vertex(m->index(i + j));

			if (vert.pos[2] < 0)
			{
				kill = true;
			}
		}
		
		if (kill)
		{
			for (size_t j = 0; j < 3; j++)
			{
				m->iPointer()[i + j] = 0;
			}

		}
	}
	
	{
		mat3x3 scale = make_mat3x3();
		mat3x3_scale(&scale, 0.6, 0.6, 1.0);
		m->setPosition(empty_vec3());
		m->rotateByMatrix(scale);
	}

	m->setColour(0.5, 0.0, 0);
	
	for (size_t i = 0; i < m->vertexCount(); i++)
	{
		Helen3D::Vertex tmp = m->vertex(i);
		vec3 v = vec_from_pos(tmp.pos);
		std::vector<Distance> closest;
		closestPoints(v, &closest);

		double wAve = 0;
		double zAve = 0;
		for (size_t j = 0; j < closest.size(); j++)
		{
			double z = closest[j].p.z;
			double d = closest[j].l;
			double weight = 1/d;
			
			wAve += weight;
			zAve += z * weight;
		}
		
		zAve /= wAve;
		m->vPointer()[i].pos[2] = zAve;
	}
	
	m->changeToTriangles();
	std::string v = Contour_vsh();
	std::string f = Contour_fsh();
	m->changeProgram(v, f);
}


void Tent::closestPoints(vec3 v, std::vector<Distance> *results)
{
	std::vector<Distance> dists;

	for (size_t i = 0; i < _ps.size(); i++)
	{
		vec3 d = _ps[i] - v;
		d.z = 0;
		double l = vec3_sqlength(d);
		Distance dist;
		dist.p = _ps[i];
		dist.l = l;
		dists.push_back(dist);
	}

	std::sort(dists.begin(), dists.end(), Tent::compare_distance);
	
	for (size_t i = 0; i < 5 && i < dists.size(); i++)
	{
		results->push_back(dists[i]);
	}
}
