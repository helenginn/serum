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

class Mutation;
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
	void writeResultVectors(std::string filename);
	void populateNames(char ****ptr);
	
	void setDimension(int dim)
	{
		_dim = dim;
		if (_dim <= 0) _dim = 1;
	}
	
	int dimensions()
	{
		return _dim;
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

	static double resultForDirection(double *dir);
	static double resultForVector(double *dir1, double *dir2);
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
	std::vector<double> _sums, _sumsqs;
	std::string _filename;

	std::vector<Strain *> _strains;
	std::vector<Serum *> _sera;
	std::vector<Mutation *> _muts;

	std::map<std::string, Strain *> _name2Strain;
	std::map<std::string, Serum *> _name2Serum;

	std::vector<Any *> _anys;

	int _count;
	bool _refineOffset;
	bool _refineStrength;
	double _scale;
	static double *_scratch;
	static int _dim;
};

#endif
