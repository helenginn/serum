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

#ifndef __serum__Strain__
#define __serum__Strain__

#include <map>

class Loader;
class Mutation;
class Serum;

class Strain
{
public:
	Strain(std::string name);
	
	void setList(std::string list, std::vector<Mutation *> &unique);
	std::vector<double> generateVector(std::vector<Mutation *> &full);
	
	double vectorCompare(Strain *str);

	const std::string &name()
	{
		return _name;
	}
	
	double strength()
	{
		return _strength;
	}
	
	double *strengthPtr()
	{
		return &_strength;
	}
	
	bool hasMutation(Mutation *mut);
	
	void clearPrecalculated()
	{
		_strainVecs.clear();
	}
	
	void setLoader(Loader *l)
	{
		_loader = l;
	}
	
	void needsRefresh(bool r = true)
	{
		_refresh = r;
	}
	
	const std::vector<double> &mutVector()
	{
		return _vec;
	}
	
	const std::vector<double> &direction()
	{
		return _dir;
	}
	
	size_t serumCount()
	{
		return _sera.size();
	}
	
	Serum *serum(int i)
	{
		return _sera[i];
	}
	
	void addSerum(Serum *s, double val, bool free = false)
	{
		_sera.push_back(s);
		_values[s] = val;
		_frees[s] = free;
	}
	
	double serumValue(Serum *s)
	{
		return _values[s];
	}
	
	bool freeSerum(Serum *s)
	{
		return _frees[s];
	}
	
	bool hasSerum(Serum *s)
	{
		return _values.count(s);
	}

	static void staticFindPosition(void *object);
	void findPosition();
private:
	std::string _name;
	std::vector<Mutation *> _list;
	std::vector<double> _vec;
	
	std::map<Strain *, std::vector<double>> _strainVecs;

	Loader *_loader;
	std::vector<Serum *> _sera;
	std::map<Serum *, double> _values;
	std::map<Serum *, bool> _frees;
	double _strength;
	bool _refresh;
	std::vector<double> _dir;
};

#endif
