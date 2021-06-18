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

#ifndef __Mutation__Mutation__
#define __Mutation__Mutation__

#include <string>
#include <vector>

class Strain;

class Mutation
{
public:
	Mutation(std::string mut, int dims = 3);

	std::string str()
	{
		return _mut;
	}
	
	bool operator==(Mutation &mut)
	{
		return mut.str() == str();
	}

	bool operator>(Mutation &mut)
	{
		return mut.str() > str();
	}
	
	void addStrain(Strain *strain)
	{
		_strains.push_back(strain);
	}
	
	double *scalarPtr(int i)
	{
		return &_vec[i];
	}
	
	double scalar(int i)
	{
		return _vec[i];
	}
	
	void randomiseVector();
	
	std::vector<double> &asVector();
	
	static void refresh(void *object);
private:
	std::string _mut;
	std::vector<Strain *> _strains;
	std::vector<double> _vector;
	int _dim;
	double *_vec;
};

#endif
