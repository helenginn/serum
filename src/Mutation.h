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

#include <hcsrc/vec3.h>
#include <hcsrc/mat3x3.h>
#include <string>
#include <vector>
#include <map>
#include <QTreeWidgetItem>

class Strain;

class Mutation : public QTreeWidgetItem
{
public:
	Mutation(std::string mut, int dims = 3);
	
	void reset();
	void resetOriginal();

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
	
	void saveVector(int idx, bool imp = false);
	void loadSaved(int i);
	
	int trials()
	{
		return _saved.size();
	}
	
	std::vector<double> &trial(int i)
	{
		return _saved[i];
	}
	
	void silence(bool s)
	{
		_silenced = s;
	}
	
	bool silenced()
	{
		return _silenced;
	}
	
	std::vector<double> &savedVector(int idx)
	{
		return _saved[idx];
	}
	
	std::vector<double> &savedOrigin(int idx)
	{
		return _savedOrigs[idx];
	}
	
	double *easePtr()
	{
		return &_ease;
	}
	
	double averageEase(double *stdev = NULL);
	
	double ease()
	{
		return _ease;
	}
	
	double *scalarPtr(int i)
	{
		return &_vec[i];
	}
	
	double scalar(int i)
	{
		return _vec[i];
	}
	
	static double origin(int i)
	{
		return _original[i];
	}
	
	double *originPtr(int i)
	{
		return &_original[i];
	}
	
	mat3x3 tensor()
	{
		return _tensor;
	}
	
	vec3 centre()
	{
		return _centre;
	}
	
	int residue()
	{
		return _res;
	}
	
	void incrementUses()
	{
		_uses++;
	}
	
	bool isUsed()
	{
		return _uses > 0;
	}

	void randomiseVector();
	
	void calculateCloud();
	
	int hash();
	
	std::vector<double> &asVector(double *from);
	double crossCorrelation(int i, int j);
	void rotateCurrent(double **rot);
	void subtract(std::vector<double> sub);
	
	static void refresh(void *object);
	static std::vector<double> centroid(std::vector<Mutation *> list);
private:
	void findWeight();
	std::string _mut;
	std::vector<Strain *> _strains;
	std::vector<double> _vector;
	std::vector<double> _eases;
	static double *_original;
	static int _dim;
	mat3x3 _tensor;
	vec3 _centre;
	double _weight;
	double *_vec;
	bool _silenced;
	double _ease;
	int _res;
	int _best;
	int _uses;
	std::map<int, std::vector<double> > _savedOrigs;
	std::map<int, std::vector<double> > _saved;
};

#endif
