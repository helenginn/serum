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

#include "Mutation.h"
#include "Strain.h"
#include "Loader.h"

Mutation::Mutation(std::string mut, int dims) : QTreeWidgetItem(0)
{
	_mut = mut;
	_silenced = false;
	_best = 0;
	_dim = dims;
	_vec = new double[dims];
	for (size_t i = 0; i < dims; i++)
	{
		_vec[i] = 0;
	}

	setText(0, QString::fromStdString(_mut));
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
		_vec[j] = ((rand() / (double)RAND_MAX) / 4);
		
		if (_silenced)
		{
			_vec[j] = 0;
		}
	}
}

void Mutation::loadSaved(int i)
{
	std::vector<double> load = savedVector(i);

	for (size_t i = 0; i < load.size(); i++)
	{
		_vec[i] = load[i];
	}
}

std::vector<double> &Mutation::asVector()
{
	_vector.resize(_dim);
	for (size_t i = 0; i < _dim; i++)
	{
		_vector[i] = fabs(_vec[i]);
	}

	return _vector;
}

void Mutation::saveVector(int idx, bool imp)
{
	if (idx < 0)
	{
		idx = _saved.size();
	}
	std::vector<double> mut = asVector();
	_saved[idx] = mut;
	
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
		return _weight * _vec[i] * _vec[j];
	}
	
	return _weight * _vec[i] * _saved[_best][j];
}

void Mutation::rotateCurrent(double **rot)
{
	asVector();
	
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
	results.resize(Loader::dimensions());
	double weights = 0;
	
	for (size_t i = 0; i < list.size(); i++)
	{
		Mutation *m = list[i];
		
		for (size_t j = 0; j < Loader::dimensions(); j++)
		{
			results[j] += m->_vec[j] * m->_weight;
			weights += m->_weight;
		}
	}

	for (size_t j = 0; j < Loader::dimensions(); j++)
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
