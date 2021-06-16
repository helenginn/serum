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

#ifndef __serum__Loader__
#define __serum__Loader__

#include <string>
#include <vector>
#include <complex>
#include <map>
#include <QObject>

class Strain;
class Serum;
class Any;

typedef std::complex<double> Complex;

class Loader : public QObject
{
Q_OBJECT
public:
	Loader();

	void load(std::string filename);
	void refine();
	void refineLoop(int count);
	void populateRaw(double ***ptr, int type = 0);
	void writeOut(std::string filename, int type);
	void antigenicity(std::string filename);
	void reorderMutations();
	void updateSums();
	void setModelToAverage();
	void removeUnused();
	void writeResultVectors(std::string filename);
	void populateNames(char ****ptr);
	
	void setDimension(int dim)
	{
		_dim = dim;
		if (_dim <= 0) _dim = 1;
	}
	
	void setRefineOffset(bool refine)
	{
		_refineOffset = refine;
	}
	
	void setRefineStrength(bool refine)
	{
		_refineStrength = refine;
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
	
	typedef enum 
	{
		ModelSimple,
		ModelVector
	} ModelType;

	static double resultForDirection(std::vector<double> dir);
	static double resultForVector(std::vector<double> dir1, 
	                              std::vector<double> dir2);
	
	void setModelType(ModelType type)
	{
		_type = type;
	}
signals:
	void resultReady();
	void update();
private:
	void prepare();
	void fillInZeros();
	void defineStrain(std::string strain);
	void defineChallenge(std::string str);
	void defineSerum(std::string str);
	void degenerateSummary();
	void prepareVectors();
	void addResult();
	void enablePairs(std::string line);
	double modelForPair(Strain *strain, Serum *serum);
	double vectorModelForPair(Strain *strain, Serum *serum);
	double score();
	double simpleScore(double weight);
	double vectorScore();
	double resultForResidue(int i);
	double restraintForResidue(int i);
	static double gradientRefresh(void *object);
	
	static double sscore(void *object)
	{
		double val =static_cast<Loader *>(object)->score();
		return val;
	}

	std::string _results;
	std::vector<double> _reals;
	std::vector<double> _sums, _sumsqs;
	std::vector<std::vector<double> > _scales, _sumScales, _sumScaleSqs;
	std::string _filename;
	std::vector<Strain *> _strains;
	std::vector<Serum *> _sera;
	std::map<std::string, Strain *> _name2Strain;
	std::map<std::string, Serum *> _name2Serum;
	std::map<Strain *, std::map<Serum *, double> > _strainMap;
	std::map<Strain *, std::map<Serum *, bool> > _freeMap;
	std::vector<std::string> _muts;
	std::vector<Any *> _anys;
	int _count;
	ModelType _type;
	bool _refineOffset;
	bool _refineStrength;
	double _scale;
	static int _dim;
};

#endif
