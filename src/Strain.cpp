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
#include "Strain.h"
#include "Loader.h"
#include "Mutation.h"
#include <iomanip>

Strain::Strain(std::string name)
{
	_refresh = true;
	_name = name;
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
			m = new Mutation(contents[i], _loader->dimensions());
			unique.push_back(m);
		}

		m->addStrain(this);
		_list.push_back(m);
	}

	std::sort(_list.begin(), _list.end(), std::greater<Mutation *>());
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
	_vec = vec;

	return vec;
}

void Strain::staticFindPosition(void *object)
{
	static_cast<Strain *>(object)->needsRefresh();
}

void Strain::findPosition()
{
	if (_refresh == false)
	{
		return;
	}

	int dim = _loader->dimensions();
	_dir = std::vector<double>(dim, 0);
	
	for (size_t i = 0; i < _list.size(); i++)
	{
		for (size_t j = 0; j < dim; j++)
		{
			_dir[j] += _list[i]->scalar(j);
		}
	}
	
	_refresh = false;
}

double Strain::vectorCompare(Strain *str) 
{
	if (str == this)
	{
		return 0;
	}
	
	return Loader::resultForVector(&_dir[0], &str->_dir[0]);
}
