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
#include <hcsrc/vec3.h>
#include <hcsrc/mat3x3.h>
#include <QTreeWidgetItem>

class Loader;
class Table;
class Mutation;
class Serum;

class Strain : public QTreeWidgetItem
{
public:
	Strain(std::string name);
	
	void setList(std::string list, std::vector<Mutation *> &unique);
	std::vector<double> generateVector(std::vector<Mutation *> &full);
	
	double vectorCompare(Strain *str);
	void sort();
	bool sameAsStrain(Strain *s);

	const std::string &name()
	{
		return _name;
	}
	
	double ease();

	void saveEase();
	
	double *easePtr()
	{
		return &_ease;
	}
	
	double strength()
	{
		return _strength;
	}
	
	double *strengthPtr()
	{
		return &_strength;
	}
	
	void setDefaultOffset(double strength)
	{
		_default = strength;
	}
	
	double offset()
	{
		return _offset;
	}
	
	double *offsetPtr()
	{
		return &_offset;
	}
	
	double *importancePtr(int j)
	{
		return &_importance[j];
	}

	bool hasMutation(Mutation *mut);
	
	void reset();
	
	void setLoader(Loader *l)
	{
		_loader = l;
	}
	
	void needsRefresh(bool r = true)
	{
		_refresh = r;
	}
	
	const std::vector<double> &direction()
	{
		return _dir;
	}
	
	double defaultOffset()
	{
		return _default;
	}
	
	size_t serumCount()
	{
		return _sera.size();
	}
	
	Serum *serum(int i)
	{
		return _sera[i];
	}
	
	void addSerum(Serum *s, Table *table, bool free = false);
	
	double serumValue(Serum *s);
	double serumScale(Serum *s);
	
	bool freeSerum(Serum *s)
	{
		return _frees[s];
	}
	
	bool hasSerum(Serum *s)
	{
		return _values.count(s);
	}
	
	size_t mutationCount()
	{
		return _list.size();
	}
	
	Mutation *mutation(int i)
	{
		return _list[i];
	}
	
	std::vector<Mutation *> &list()
	{
		return _list;
	}

	static void staticFindPosition(void *object);
	void findPosition();
	vec3 aveDirection();
	void calculateCloud();

	std::string summary();
	
	int trials();
	
	double averageEase();
	
	bool isFree()
	{
		return _free;
	}
	
	void setFree(bool f)
	{
		_free = f;
	}
	
	mat3x3 &tensor()
	{
		return _tensor;
	}
private:
	std::string _name;
	std::vector<Mutation *> _list;
	
	std::map<Strain *, std::vector<double>> _strainVecs;

	Loader *_loader;
	std::vector<Serum *> _sera;
	std::map<Serum *, std::vector<Table *> > _values;
	std::map<Serum *, bool> _frees;
	std::vector<double> _eases;
	std::vector<double> _selfEases;
	double _default;
	double _strength;
	double _offset;
	double _ease;
	bool _refresh;
	bool _free;
	mat3x3 _tensor;
	std::vector<double> _dir, _scaledDir;
	std::vector<double> _importance;
};

#endif
