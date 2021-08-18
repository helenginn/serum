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
#include "Table.h"
#include "Serum.h"
#include "Strain.h"

Table::Table(std::string name)
{
	_name = name;
	_count = 0;
	_scale = 1;
}

void Table::addValue(Strain *strain, Serum *serum, double value)
{
	if (std::find(_strains.begin(), _strains.end(), strain) == _strains.end())
	{
		_strains.push_back(strain);
	}

	if (std::find(_sera.begin(), _sera.end(), serum) == _sera.end())
	{
		_sera.push_back(serum);
	}
	
	_values[strain][serum].push_back(value);
	
	int num = _values[strain][serum].size();
	double ave = 0;
	for (size_t i = 0; i < num; i++)
	{
		ave += _values[strain][serum][i];
	}
	ave /= num;
	_averages[strain][serum] = ave;
	_count++;

	strain->addSerum(serum, this);
	
	Challenge ch;
	ch.strain = strain;
	ch.serum = serum;
	ch.value = value;
	_challenges.push_back(ch);
}

double Table::valueFor(Strain *strain, Serum *serum)
{
	if (_averages.count(strain) && _averages[strain].count(serum))
	{
		return _averages[strain][serum];
	}
	
	return NAN;
}

void Table::reset()
{
	_scale = 1;
}

double Table::challenge(int i, Strain **strainPtr, Serum **serumPtr)
{
	Challenge &ch = _challenges[i];
	*strainPtr = ch.strain;
	*serumPtr = ch.serum;

	return ch.value;
}
