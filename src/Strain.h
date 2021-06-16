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

class Strain
{
public:
	Strain(std::string name);
	
	void setList(std::string list);
	void addToMuts(std::vector<std::string> &unique);
	std::vector<double> generateVector(std::vector<std::string> &full);
	double compareToStrain(Strain *str, std::vector<double> &scales,
	                        std::vector<std::string> &coops);

	
	double vectorCompare(Strain *str, std::vector<std::vector<double> > &scales);

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
	
	bool hasMutation(std::string mut);
	
	void clearPrecalculated()
	{
		_strainVecs.clear();
	}
	
	const std::vector<double> &mutVector()
	{
		return _vec;
	}
	
	const std::vector<double> &direction()
	{
		return _dir;
	}

	static bool hasCooperative(Strain *a, Strain *b, std::string coop);

	std::vector<std::string> combinedLists(Strain *other);
	void findPosition(std::vector<std::vector<double> > &scales);
private:
	std::vector<double> scaleVector(std::vector<double> &scales,
	                                std::vector<std::string> &coops,
	                                Strain *other);
	std::string _name;
	std::vector<std::string> _list;
	std::vector<double> _vec;
	
	std::map<Strain *, std::vector<double>> _strainVecs;

	double _strength;
	std::vector<double> _dir;
};

#endif
