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
#include <iomanip>

Strain::Strain(std::string name)
{
	_name = name;
}

void Strain::setList(std::string list)
{
	_list = split(list, ',');

	for (size_t i = 0; i < _list.size(); i++)
	{
		trim(_list[i]);
	}

	std::sort(_list.begin(), _list.end(), std::greater<std::string>());
}

bool Strain::hasCooperative(Strain *a, Strain *b, std::string coop)
{
	if (a == NULL || b == NULL || coop.find('+') == std::string::npos)
	{
		return false;
	}

	std::vector<std::string> each = split(coop, '+');
	
	if (each.size() != 2)
	{
		return false;
	}

	std::vector<std::string>::iterator it_a0, it_a1, it_b0, it_b1;

	it_a0 = std::find(a->_list.begin(), a->_list.end(), each[0]);
	it_b0 = std::find(b->_list.begin(), b->_list.end(), each[0]);
	it_a1 = std::find(a->_list.begin(), a->_list.end(), each[1]);
	it_b1 = std::find(b->_list.begin(), b->_list.end(), each[1]);

	if (it_a0 != a->_list.end() && it_a1 != a->_list.end())
	{
		return false;
	}
	else if (it_b0 != b->_list.end() && it_b1 != b->_list.end())
	{
		return false;
	}
	else if ((it_a0 != a->_list.end() && it_b1 != b->_list.end()) ||
	    (it_a1 != a->_list.end() && it_b0 != b->_list.end()))
	{
		return true;
	}
	
	return false;
}

void Strain::addToMuts(std::vector<std::string> &unique)
{
	std::cout << "Strain " << _name << " with " << std::flush;
	std::vector<std::string>::iterator it;

	for (size_t i = 0; i < _list.size(); i++)
	{
		std::cout << _list[i] << " " << std::flush;
		it = std::find(unique.begin(), unique.end(), _list[i]);

		if (it == unique.end())
		{
			unique.push_back(_list[i]);
		}
	}

	if (_list.size() == 0)
	{
		std::cout << "[N/A]";
	}

	std::cout << std::endl;
}

bool Strain::hasMutation(std::string mut)
{
	return (std::find(_list.begin(), _list.end(), mut) != _list.end());
}

std::vector<double> Strain::generateVector(std::vector<std::string> &full)
{
	std::vector<double> vec = std::vector<double>(full.size(), 0);
	std::cout << std::setw(20) << _name << ": ";

	std::vector<std::string>::iterator it;
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

std::vector<std::string> Strain::combinedLists(Strain *other)
{
	std::vector<std::string> all = _list;
	all.reserve(_list.size() + other->_list.size());
	all.insert(all.end(), other->_list.begin(), other->_list.end());

	return all;
}

std::vector<double> Strain::scaleVector(std::vector<double> &scales,
                                        std::vector<std::string> &muts,
                                        Strain *other)
{
	std::vector<double> results;
	
	if (_strainVecs.count(other) == 0)
	{
		for (size_t i = 0; i < _vec.size(); i++)
		{
			double c = _vec[i];

			if (c <= 0 && hasCooperative(this, other, muts[i]))
			{
				c = 1;
			}

			results.push_back(c);
		}

		_strainVecs[other] = results;
	}
	else
	{
		results = _strainVecs[other];
	}

	for (size_t i = 0; i < scales.size(); i++)
	{
		results[i] *= scales[i];
	}

	return results;
}

double Strain::compareToStrain(Strain *str, std::vector<double> &scales,
                                std::vector<std::string> &muts)
{
	std::vector<double> mine = scaleVector(scales, muts, str);
	std::vector<double> other = str->scaleVector(scales, muts, str);

	double sum = 0;
	for (size_t i = 0; i < scales.size(); i++)
	{
		bool pos = scales[i] > 0;
		double mult = (mine[i] - other[i]);
		mult *= mult;
		if (!pos)
		{
			mult *= -1;
		}
		sum += mult;
	}
	
	return sum;
}

void Strain::findPosition(std::vector<std::vector<double> > &scales)
{
	if (scales.size() == 0) return;
	std::vector<double> me = mutVector();
	
	int dim = scales[0].size();
	_dir = std::vector<double>(dim, 0);
	
	for (size_t i = 0; i < scales.size(); i++)
	{
		std::vector<double> mymut = scales[i];

		for (size_t j = 0; j < dim; j++)
		{
			_dir[j] += mymut[j] * me[i];
		}
	}
}

double Strain::vectorCompare(Strain *str, 
                           std::vector<std::vector<double> > &scales)
{
	if (str == this)
	{
		return 0;
	}
	
	return Loader::resultForVector(_dir, str->_dir);
}
