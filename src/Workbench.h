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

#ifndef __serum__Workbench__
#define __serum__Workbench__

#include <string>
#include <vector>
#include <complex>
#include <map>
#include <QObject>

class Mutation;
class Strain;
class Serum;
class Table;
class Any;

class Workbench : public QObject
{
Q_OBJECT
public:
	Workbench();

	void refine(int cycle);
	void populateRaw(double ***ptr, int type = 0);
	void writeOut(std::string filename, int type);
	void antigenicity(std::string filename);
	void reorderMutations();
	void updateSums();
	void setModelToAverage();
	void writeResultVectors(std::string filename);
	void populateNames(char ****ptr);
	
	void addMutations(std::vector<Mutation *> mutations);
	void addStrains(std::vector<Strain *> strains);
	void addSera(std::vector<Serum *> sera);
	void addTables(std::vector<Table *> tables);

	void makeStrainFree(std::string str);
	
	const std::vector<Mutation *> &mutations()
	{
		return _muts;
	}

	const std::vector<Strain *> &strains()
	{
		return _strains;
	}
	
	Strain *strain(int i)
	{
		return _strains[i];
	}
	
	Strain *strain(std::string name)
	{
		if (_name2Strain.count(name) == 0)
		{
			return NULL;
		}

		return _name2Strain[name];
	}
	
	void setDimension(int dim)
	{
		_dim = dim;
		if (_dim < 0) _dim = 1;
	}
	
	static int dimensions()
	{
		return _dim;
	}
	
	void setRefineEase(bool ease)
	{
		_refineEase = ease;
	}
	
	void setRefineOffset(bool refine)
	{
		_refineOffset = refine;
	}
	
	void setRefineStrainStrength(bool refine)
	{
		_refineStrainStrength = refine;
	}
	
	void setRefineStrainOffset(bool refine)
	{
		_refineStrainOffset = refine;
	}

	void setRefineStrength(bool refine)
	{
		_refineStrength = refine;
	}

	void setRefineImportance(bool refine)
	{
		_refineImportance = refine;
	}
	
	void perResidueAntigenicity(std::string filename);
	
	void setScale(double scale)
	{
		_scale = scale;
	}
	
	size_t strainCount()
	{
		return _strains.size();
	}
	
	size_t serumCount()
	{
		return _sera.size();
	}
	
	void displaySettings(std::string filename);

	static double resultForDirection(double *dir);
	static double resultForVector(double *dir1, double *dir2);
	void markStrains(std::string str);
	void markSera(std::string str);
	void setup();

	void clusterMutations(std::string filename);
	double modelForPair(Strain *strain, Serum *serum);
	double aveModel(Strain *strain, Serum *serum);
	Strain *nextStrain(double *score);
	std::vector<Strain *> getStrains(std::string list);
public slots:
	void refineLoop();
signals:
	void resultReady();

private:
	void startOffsets();
	void rotations();
	void bestRotation();
	void prepare();
	void prepareVectors();
	void addResult();
	double vectorModelForPair(Strain *strain, Serum *serum);
	double score();
	double simpleScore(double weight);
	double vectorScore();
	double resultForResidue(int i);
	double restraintForResidue(int i);
	void degenerateSummary();
	static double gradientRefresh(void *object);
	
	static double sscore(void *object)
	{
		double val =static_cast<Workbench *>(object)->score();
		return val;
	}

	std::string _results;
	std::vector<double> _sums, _sumsqs;

	std::vector<Strain *> _strains;
	std::vector<Serum *> _sera;
	std::vector<Mutation *> _muts;
	std::vector<Table *> _tables;

	std::map<std::string, Strain *> _name2Strain;
	std::map<std::string, Serum *> _name2Serum;

	std::vector<Any *> _anys;

	int _count;
	bool _shouldSetup;
	bool _refineEase;
	bool _refineOffset;
	bool _refineStrength;
	bool _refineImportance;
	bool _refineStrainStrength;
	bool _refineStrainOffset;
	double _scale;
	double _best;
	static double *_scratch;
	std::vector<double> _scores;
	static int _dim;
};

#endif
