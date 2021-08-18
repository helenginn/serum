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

#include <hcsrc/FileReader.h>
#include <libsrc/Anisotropicator.h>
#include "Table.h"
#include "Strain.h"
#include "Serum.h"
#include "Workbench.h"
#include "Mutation.h"
#include "Loader.h"
#include <iomanip>

Strain::Strain(std::string name)
{
	_free = false;
	_refresh = true;
	_name = name;
	setText(0, QString::fromStdString(name));
	_strength = 1;
	_default = 0;
	_offset = 0;
	_ease = 0;
	int dim = Workbench::dimensions();
	_importance = std::vector<double>(dim, 1);
}

void Strain::saveEase()
{
	_eases.push_back(ease());
	_selfEases.push_back(_ease);
}

void Strain::reset()
{
	int dim = Workbench::dimensions();
	_strength = 1;
	_ease = 0;
	_offset = _default;
	_importance = std::vector<double>(dim, 1);
	needsRefresh();
}

void Strain::setList(std::string list, std::vector<Mutation *> &unique)
{
	_list.clear();

	std::vector<std::string> contents = split(list, ',');

	for (size_t i = 0; i < contents.size(); i++)
	{
		trim(contents[i]);
		
		Mutation *m = NULL;
		
		for (size_t j = 0; j < unique.size(); j++)
		{
			if (unique[j]->str() == contents[i])
			{
				m = unique[j];
			}
		}
		
		if (m == NULL)
		{
			m = new Mutation(contents[i], Workbench::dimensions());
			
			if (!_loader->isAcceptable(m->residue()))
			{
				continue;
			}

			unique.push_back(m);
		}

		m->addStrain(this);
		_list.push_back(m);
	}

	sort();
}

bool Strain::hasMutation(Mutation *mut)
{
	return (std::find(_list.begin(), _list.end(), mut) != _list.end());
}

std::vector<double> Strain::generateVector(std::vector<Mutation *> &full)
{
	std::vector<double> vec = std::vector<double>(full.size(), 0);
	std::cout << std::setw(20) << _name << ": ";

	std::vector<Mutation *>::iterator it;
	for (size_t i = 0; i < full.size(); i++)
	{
		it = std::find(_list.begin(), _list.end(), full[i]);
		
		if (it != _list.end())
		{
			vec[i] = 1;
		}
		
		std::cout << vec[i] << " " << std::flush;
	}
	
	std::cout << std::endl;

	return vec;
}

void Strain::staticFindPosition(void *object)
{
	static_cast<Strain *>(object)->needsRefresh();
}

vec3 Strain::aveDirection()
{
	vec3 ave = empty_vec3();

	for (size_t i = 0; i < _list.size(); i++)
	{
		_list[i]->calculateCloud();
		ave += _list[i]->centre();
	}

	return ave;
}

void Strain::findPosition()
{
	if (_refresh == false)
	{
		return;
	}

	int dim = Workbench::dimensions();
	_dir = std::vector<double>(dim, 0);
	_scaledDir = std::vector<double>(dim, 0);

	for (size_t j = 0; j < dim; j++)
	{
		_dir[j] = Mutation::origin(j);
	}
	
	for (size_t i = 0; i < _list.size(); i++)
	{
		for (size_t j = 0; j < dim; j++)
		{
			_dir[j] += _list[i]->scalar(j);
		}
	}
	
	for (size_t i = 0; i < dim; i++)
	{
		_scaledDir[i] = _importance[i] * _dir[i];
	}
	
	_refresh = false;
}

double Strain::vectorCompare(Strain *str) 
{
	if (str == this)
	{
		return 0;
	}
	
	findPosition();
	str->findPosition();
	
	double val = Workbench::resultForVector(&_scaledDir[0], &str->_dir[0]);
	val -= str->ease();
	return val;
}

void Strain::addSerum(Serum *s, Table *t, bool free)
{
	s->addStrain(this);
	_sera.push_back(s);
	_values[s].push_back(t);
	_frees[s] = free;
}

double Strain::serumScale(Serum *s)
{
	if (_values.count(s) && _values[s].size() > 0)
	{
		Table *t = _values[s][0];
		return t->scale();
	}
	
	return NAN;
}

double Strain::serumValue(Serum *s)
{
	if (_values.count(s) && _values[s].size() > 0)
	{
		Table *t = _values[s][0];
		return t->valueFor(this, s);
	}
	
	return NAN;
}

std::string Strain::summary()
{
	std::string sum;
	for (size_t i = 0; i < _list.size(); i++)
	{
		sum += _list[i]->str() + ", ";
	}
	
	if (_list.size() > 0)
	{
		sum.pop_back();
		sum.pop_back();
	}

	return sum;
}

void Strain::sort()
{
	std::sort(_list.begin(), _list.end(), std::greater<Mutation *>());
}

bool Strain::sameAsStrain(Strain *s)
{
	sort();
	s->sort();
	
	return (_list == s->_list);
}

void Strain::calculateCloud()
{
	if (_list.size() == 0)
	{
		return;
	}

	int dim = Workbench::dimensions();
	std::vector<std::vector<double> > all;

	for (size_t j = 0; j < _list[0]->trials(); j++)
	{
		std::vector<double> dir = std::vector<double>(dim, 0);

		for (size_t i = 0; i < _list.size(); i++)
		{
			std::vector<double> vals = _list[i]->trial(j);
			for (size_t k = 0; k < dim; k++)
			{
				dir[k] += vals[k];
			}
		}
		
		all.push_back(dir);
	}

	vec3 centre = aveDirection();

	std::vector<vec3> points;
	for (size_t i = 0; i < all.size(); i++)
	{
		vec3 point = empty_vec3();
		for (size_t j = 0; j < dim && j < 3; j++)
		{
			double val = all[i][j];
			double mid = *(&centre.x + j);
			*(&point.x + j) = val - mid;
		}

		points.push_back(point);
	}

	Anisotropicator aniso;
	aniso.setPoints(points);
	_tensor = aniso.basis();

	for (int i = 0; i < 3; i++)
	{
		vec3 ax = mat3x3_axis(_tensor, i);
		double l = vec3_length(ax);
		double scale = sqrt(l) / sqrt(trials());
		if (scale == scale)
		{
			vec3_set_length(&ax, scale);
			mat3x3_set_axis(&_tensor, i, ax);
		}
	}
	
	for (size_t i = 0; i < 9; i++)
	{
		if (_tensor.vals[i] != _tensor.vals[i])
		{
			_tensor.vals[i] = 0;
		}
	}
}

int Strain::trials()
{
	if (_list.size() == 0)
	{
		return 0;
	}
	
	return _list[0]->trials();
}

double Strain::averageEase()
{
	double sum = 0; 
	for (size_t i = 0; i < _eases.size(); i++)
	{
		sum += _eases[i];
	}
	
	sum /= (double)_eases.size();
	return sum;
}

double Strain::ease()
{
	double sum = 0; 
	for (size_t i = 0; i < _list.size(); i++)
	{
		sum += _list[i]->ease();
	}
	
	sum += _ease;
	
	return sum;
}
