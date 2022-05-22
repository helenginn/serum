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

#include <iostream>
#include <libsrc/Anisotropicator.h>
#include "Mutation.h"
#include "Strain.h"
#include "Workbench.h"

double *Mutation::_original = NULL;
int Mutation::_dim = 1;

Mutation::Mutation(std::string mut, int dims) : QTreeWidgetItem(0)
{
	_mut = mut;
	_silenced = false;
	_best = 0;
	_uses = 0;
	_ease = 0;
	_dim = dims;
	
	if (_original == NULL)
	{
		_original = new double[dims];
		for (size_t i = 0; i < dims; i++)
		{
			_original[i] = 0;
		}
	}

	_vec = new double[dims];
	for (size_t i = 0; i < dims; i++)
	{
		_vec[i] = 0;
	}

	setText(0, QString::fromStdString(_mut));
	_res = -1;

	if (_mut.length() == 0)
	{
		return;
	}

	char *old = &_mut[0];
	char *pos = old;
	int count = 0;
	int val = 0;

	while (pos == old && count <= _mut.length())
	{
		count++;
		old++;
	}
	
	if (count <= _mut.length())
	{
		_res = val;
	}
}

void Mutation::resetOriginal()
{
	for (size_t i = 0; i < _dim; i++)
	{
		_original[i] = 0;
	}
}

void Mutation::reset()
{
	_ease = 0;
}

void Mutation::refresh(void *object)
{
	Mutation *mut = static_cast<Mutation *>(object);
	
	for (size_t i = 0; i < mut->_strains.size(); i++)
	{
		Strain *str = mut->_strains[i];
		
		str->needsRefresh();
	}
}

void Mutation::randomiseVector()
{
	for (size_t j = 0; j < _dim; j++)
	{
		_vec[j] = ((rand() / (double)RAND_MAX) - 0.5) / 10;
		
		if (_silenced)
		{
			_vec[j] = 0;
		}
		
		_ease = 0.05 * (rand() / (double)RAND_MAX - 0.5);
	}
}

void Mutation::loadSaved(int i)
{
	std::vector<double> load = savedVector(i);

	for (size_t i = 0; i < load.size(); i++)
	{
		_vec[i] = load[i];
	}

	load = savedOrigin(i);
	
	for (size_t i = 0; i < load.size(); i++)
	{
		_original[i] = load[i];
	}
}

std::vector<double> &Mutation::asVector(double *from)
{
	_vector.resize(_dim);
	for (size_t i = 0; i < _dim; i++)
	{
		_vector[i] = from[i];
	}

	return _vector;
}

void Mutation::saveVector(int idx, bool imp)
{
	if (idx < 0)
	{
		idx = _saved.size();
	}
	std::vector<double> mut = asVector(_vec);
	_saved[idx] = mut;
	_eases.push_back(_ease);
	std::vector<double> orig = asVector(_original);
	_savedOrigs[idx] = orig;
	
	findWeight();
	if (imp)
	{
		_best = idx;
	}
}

void Mutation::findWeight()
{
	std::vector<double> all;
	all.resize(_dim);

	for (size_t i = 0; i < trials(); i++)
	{
		std::vector<double> v = savedVector(i);
		for (size_t j = 0; j < _dim; j++)
		{
			all[j] += v[j];
		}
	}

	_weight = 0;
	for (size_t j = 0; j < _dim; j++)
	{
		_weight += all[j] * all[j];
	}

	_weight = sqrt(_weight);
}

double Mutation::crossCorrelation(int i, int j)
{
	if (_saved.count(_best) == 0)
	{
		return 1 * _vec[i] * _vec[j];
	}
	
	return 1 * _vec[i] * _saved[_best][j];
}

void Mutation::rotateCurrent(double **rot)
{
	asVector(_vec);
	
	for (size_t i = 0; i < _vector.size(); i++)
	{
		_vector[i] = 0;
		for (size_t j = 0; j < _vector.size(); j++)
		{
			_vector[i] += _vec[j] * rot[i][j];
		}
	}

	for (size_t i = 0; i < _vector.size(); i++)
	{
		_vec[i] = _vector[i];
	}
}

std::vector<double> Mutation::centroid(std::vector<Mutation *> list)
{
	std::vector<double> results;
	results.resize(Workbench::dimensions());
	double weights = 0;
	
	for (size_t i = 0; i < list.size(); i++)
	{
		Mutation *m = list[i];
		
		for (size_t j = 0; j < Workbench::dimensions(); j++)
		{
			results[j] += m->_vec[j] * m->_weight;
			weights += m->_weight;
		}
	}

	for (size_t j = 0; j < Workbench::dimensions(); j++)
	{
		results[j] /= weights;
	}
	
	return results;
}

void Mutation::subtract(std::vector<double> sub)
{
	for (size_t i = 0; i < sub.size(); i++)
	{
		_vec[i] -= sub[i];
	}
}

void Mutation::calculateCloud()
{
	std::vector<double> all;
	all.resize(_dim);

	for (size_t i = 0; i < trials(); i++)
	{
		std::vector<double> v = savedVector(i);
		for (size_t j = 0; j < _dim; j++)
		{
			all[j] += v[j];
		}
	}

	std::vector<vec3> points;
	for (size_t i = 0; i < trials(); i++)
	{
		std::vector<double> v = savedVector(i);
		v.resize(3);
		vec3 point = empty_vec3();
		for (size_t j = 0; j < _dim && j < 3; j++)
		{
			*(&point.x + j) = v[j] - all[j] / (double)trials();
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
		vec3_set_length(&ax, sqrt(l));
		mat3x3_set_axis(&_tensor, i, ax);
	}
	
	for (size_t i = 0; i < 9; i++)
	{
		if (_tensor.vals[i] != _tensor.vals[i])
		{
			_tensor.vals[i] = 0;
		}
	}

	vec3 centre = empty_vec3();
	all.resize(3);
	centre.x = all[0] / (double)trials();
	centre.y = all[1] / (double)trials();
	centre.z = all[2] / (double)trials();
	
	_centre = centre;
}

double Mutation::averageEase(double *stdev)
{
	double sum = 0; 
	double sumsq = 0;
	for (size_t i = 0; i < _eases.size(); i++)
	{
		sum += _eases[i];
		sumsq += _eases[i] * _eases[i];
	}
	
	sum /= (double)_eases.size();
	sumsq /= (double)_eases.size();
	
	if (stdev)
	{
		*stdev = sqrt(sumsq - sum * sum);
	}

	return sum;
}

int Mutation::hash()
{
	int num = 0;
	for (size_t j = 0; j < _mut.length(); j++)
	{
		num += (int)(_mut[j]) - (int)'A';
	}

	return 6 * (num % 60);
}
