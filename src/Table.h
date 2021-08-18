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

#ifndef __serum__Table__
#define __serum__Table__

#include <map>
#include <vector>

class Strain;
class Serum;

class Table
{
public:
	Table(std::string name);

	void addValue(Strain *strain, Serum *serum, double value);
	double valueFor(Strain *strain, Serum *serum);
	void reset();
	
	double scale()
	{
		return _scale;
	}
	
	double *scalePtr()
	{
		return &_scale;
	}
	
	size_t valueCount()
	{
		return _count;
	}

	size_t strainCount()
	{
		return _strains.size();
	}
	
	size_t serumCount()
	{
		return _sera.size();
	}
	
	size_t challengeCount()
	{
		return _challenges.size();
	}
	
	double challenge(int i, Strain **strainPtr, Serum **serumPtr);
	
	Serum *serum(int i)
	{
		return _sera[i];
	}
	
	Strain *strain(int i)
	{
		return _strains[i];
	}
private:
	std::vector<Strain *> _strains;
	std::vector<Serum *> _sera;
	std::map<Strain *, std::map<Serum *, std::vector<double> > > _values;
	std::map<Strain *, std::map<Serum *, double > > _averages;
	
	typedef struct 
	{
		Strain *strain;
		Serum *serum;
		double value;
	} Challenge;

	std::vector<Challenge> _challenges;

	std::string _name;
	double _scale;
	int _count;
};

#endif
