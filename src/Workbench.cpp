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

#include "Workbench.h"
#include "Settings.h"
#include "Strain.h"
#include "Mutation.h"
#include "Serum.h"
#include "Table.h"

#include <libsrc/shared_ptrs.h>
#include <hcsrc/RefinementNelderMead.h>
#include <hcsrc/RefinementLBFGS.h>
#include <hcsrc/maths.h>
#include <hcsrc/Any.h>
#include <libica/svdcmp.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <hcsrc/FileReader.h>

int Workbench::_dim = 2;
double *Workbench::_scratch = NULL;

Workbench::Workbench()
{
	_best = 1000;
	_shouldSetup = false;
	_refineEase = false;
	_refineOffset = false;
	_refineStrength = false;
	_refineImportance = false;
	_refineStrainOffset = false;
	_refineStrainStrength = false;
	_count = 0;
	_scale = 0.3;
}

void Workbench::startOffsets()
{
	std::map<Strain *, double> scores;
	std::map<Strain *, double> counts;

	for (size_t j = 0; j < serumCount(); j++)
	{
		Serum *ser = _sera[j];
		Strain *str = ser->strain();
		
		if (!str->hasSerum(ser))
		{
			continue;
		}

		double val = str->serumValue(ser);
		scores[str] += val;
		counts[str]++;
	}
	
	std::cout << scores.size() << " strains found" << std::endl;
	std::vector<Strain *> unfinished;

	for (size_t j = 0; j < strainCount(); j++)
	{
		Strain *str = _strains[j];
		if (scores.count(str) == 0)
		{
			unfinished.push_back(str);
			continue;
		}
		
		double ave = scores[str] / counts[str];
		str->setDefaultOffset(ave);
		std::cout << str->name() << " " << ave << std::endl;
	}

	for (size_t i = 0; i < unfinished.size(); i++)
	{
		Strain *next = unfinished[i];
		double sum = 0;
		double count = 0;
		
		for (size_t j = 0; j < serumCount(); j++)
		{
			if (_sera[j]->strain() != next)
			{
				continue;
			}
			
			double max = -FLT_MAX;
			for (size_t k = 0; k < strainCount(); k++)
			{
				if (!_strains[k]->hasSerum(_sera[j]))
				{
					continue;
				}

				double val = _strains[k]->serumValue(_sera[j]);
				if (val != val)
				{
					continue;
				}
				
				if (val > max)
				{
					max = val;
				}
			}

			if (max > -FLT_MAX)
			{
				count++;
				sum += max;
			}
		}
		
		sum /= count;
		
		if (sum == sum)
		{
			next->setDefaultOffset(sum);
		}
	}
}

void Workbench::prepareVectors()
{
	std::cout << "Vectors: " << std::endl;
	
	for (size_t i = 0; i < _strains.size(); i++)
	{
		std::vector<double> vecs = _strains[i]->generateVector(_muts);
	}

	std::cout << std::endl;
	
	std::cout << "All " << _muts.size() << " mutations over " << _strains.size() 
	<< " strains: " << std::endl;

	int usable = 0;
	for (size_t i = 0; i < _muts.size(); i++)
	{
		std::cout << _muts[i]->str() << " " << std::flush;
		if (!_muts[i]->silenced() && _muts[i]->isUsed())
		{
			usable++;
		}
	}

	std::cout << std::endl;
	std::cout << "Of those: " << usable << " for refinement." << std::endl;

	size_t count = 0;
	
	for (size_t i = 0; i < _tables.size(); i++)
	{
		count += _tables[i]->valueCount();
	}

	std::cout << "Total number of observations: " << count << std::endl;
}

void Workbench::degenerateSummary()
{
	std::map<std::string, std::vector<Mutation *>> map;
	for (size_t i = 0; i < _muts.size(); i++)
	{
		Mutation *m = _muts[i];
		std::string v;
		
		for (size_t j = 0; j < _strains.size(); j++)
		{
			bool has = (_strains[j]->hasMutation(m));
			v += (has ? "1" : "0");
		}

		map[v].push_back(m);
	}
	
	std::map<std::string, std::vector<Mutation *>>::iterator it;

	std::cout << std::endl;
	std::cout << "Mutations grouped by degeneracy:" << std::endl;
	for (it = map.begin(); it != map.end(); it++)
	{
		std::vector<Mutation *> v = it->second;
		
		for (size_t j = 0; j < v.size(); j++)
		{
			std::cout << v[j]->str() << ", ";
		}

		std::cout << std::endl;
	}
	
	std::cout << std::endl;
}

typedef struct
{
	Mutation *mut;
	int res;
} StrRes;

bool higher_seq(StrRes &a, StrRes &b)
{
	if (a.res == b.res)
	{
		return a.mut->str() < b.mut->str();
	}
	return (a.res < b.res);
}

void Workbench::reorderMutations()
{
	std::cout << "Reordering mutations!" << std::endl;
	std::vector<StrRes> residues;
	for (size_t i = 0; i < _muts.size(); i++)
	{
		StrRes res;
		res.mut = _muts[i];
		std::string mut = _muts[i]->str();
		char *old = &mut[0];
		char *pos = old;
		int val = 0;
		bool found = false;
		
		while ((*old < '0' || *old > '9') && *old != '\0')
		{
			old++;
		}

		if (*old != '\0')
		{
			val = strtol(old, &pos, 10);
		}
			
		res.res = val;

		residues.push_back(res);
	}

	std::sort(residues.begin(), residues.end(), higher_seq);
	
	for (size_t i = 0; i < _muts.size(); i++)
	{
		_muts[i] = residues[i].mut;
	}
}

void Workbench::refineLoop()
{
	int count = 20;
	std::cout << "*** REFINE ***" << std::endl;

	prepare();

	for (size_t i = 0; i < count; i++)
	{
		refine(i);
	}
	
	updateSums();
	
	addResult();

	emit resultReady();
}

void Workbench::addResult()
{
	std::string add = "";
	add += "result_" + f_to_str(score(), 3) + "_";
	add += i_to_str(_count) + ",";
	
	for (size_t i = 0; i < _muts.size(); i++)
	{
		double result = resultForResidue(i);
		add += f_to_str(result, 3) + ",";
	}

	add += "\n";

	_results += add;
}

void Workbench::writeResultVectors(std::string filename)
{
	std::ofstream file;
	file.open(filename);
	
	file << "result_num,";
	for (size_t i = 0; i < _muts.size(); i++)
	{
		file << _muts[i]->str() << ",";
	}
	file << std::endl;
	
	file << _results;
	file.close();
}

void Workbench::prepare()
{
	srand(_count + 1);
	std::cout << "Random number: " << rand() << std::endl;

	for (size_t i = 0; i < _strains.size(); i++)
	{
		_strains[i]->reset();
	}

	for (size_t i = 0; i < _muts.size(); i++)
	{
		_muts[i]->randomiseVector();
	}

	for (size_t i = 0; i < _sera.size(); i++)
	{
		(*_sera[i]->strengthPtr()) = 1;
		(*_sera[i]->offsetPtr()) = 0;
	}
	
	for (size_t i = 0; i < _tables.size(); i++)
	{
		_tables[i]->reset();
	}

	std::cout << "Starting score: " << score() << std::endl;
}

void Workbench::refine(int cycle)
{
	if (_shouldSetup)
	{
		setup();
	}

	score();
	RefinementLBFGSPtr n = RefinementLBFGSPtr(new RefinementLBFGS());
	n->setCycles(1000);
	
	for (size_t i = 0; i < _anys.size(); i++)
	{
		delete _anys[i];
	}
	_anys.clear();

	for (size_t i = 0; i < _muts.size(); i++)
	{
		if (_muts[i]->silenced() || !_muts[i]->isUsed())
		{
			continue;
		}

		for (size_t j = 0; j < _dim; j++)
		{
			Any *any_real = new Any(_muts[i]->scalarPtr(j), 0.5);
			_anys.push_back(any_real);
			any_real->setRefresh(Mutation::refresh, _muts[i]);

			std::string mut = _muts[i]->str();

			n->addParameter(any_real, Any::get, Any::set, 0.2, 0.01, 
			                mut + "_" + i_to_str(j));
		}

		if (_refineEase && cycle > 10)
		{
			Any *any_real = new Any(_muts[i]->easePtr(), 0.5);
			_anys.push_back(any_real);

			std::string str = _muts[i]->str();

			n->addParameter(any_real, Any::get, Any::set, 0.2, 0.01, 
			                str + "_ease");
		}
	}

	for (size_t i = 0; i < _tables.size() && cycle > 10 && false; i++)
	{
		Any *any_real = new Any(_tables[i]->scalePtr());
		_anys.push_back(any_real);

		std::string str = "table_" + i_to_str(i);
		n->addParameter(any_real, Any::get, Any::set, 0.5, 0.01, str);
	}

	for (size_t i = 0; i < _strains.size(); i++)
	{
		if (_strains[i]->childrenCount() == 0)
		{
			continue;
		}

		for (size_t j = 0; j < _dim && _refineImportance; j++)
		{
			Any *any_real = new Any(_strains[i]->importancePtr(j));
			_anys.push_back(any_real);
			any_real->setRefresh(Strain::staticFindPosition, _strains[i]);

			std::string str = _strains[i]->name();

			n->addParameter(any_real, Any::get, Any::set, 0.2, 0.02, 
			                str + "_imp_" + i_to_str(j));
		}
		
		if (_refineStrainOffset)
		{
			Any *any_real = new Any(_strains[i]->offsetPtr());
			_anys.push_back(any_real);
			any_real->setRefresh(Strain::staticFindPosition, _strains[i]);

			std::string str = _strains[i]->name();

			n->addParameter(any_real, Any::get, Any::set, 0.2, 0.01, 
			                str + "_offset");
		}
		
		if (_refineStrainStrength)
		{
			Any *any_real = new Any(_strains[i]->strengthPtr());
			_anys.push_back(any_real);
			any_real->setRefresh(Strain::staticFindPosition, _strains[i]);

			std::string str = _strains[i]->name();

			n->addParameter(any_real, Any::get, Any::set, 0.2, 0.01, 
			                str + "_strength");
		}
	}

	for (size_t i = 0; i < _sera.size(); i++)
	{
		if (_sera[i]->strainCount() == 0)
		{
			continue;
		}

		std::string name = _sera[i]->name();

		if (_refineStrength)
		{
			double *strength = _sera[i]->strengthPtr();
			Any *any = new Any(strength);
			_anys.push_back(any);
			any->setRefresh(Strain::staticFindPosition, _sera[i]->strain());

			n->addParameter(any, Any::get, Any::set, 0.5, 0.01, name + "_s");
		}

		if (_refineOffset)
		{
			double *offset = _sera[i]->offsetPtr();
			Any *any = new Any(offset);
			_anys.push_back(any);
			any->setRefresh(Strain::staticFindPosition, _sera[i]->strain());

			n->addParameter(any, Any::get, Any::set, 0.2, 0.01, name + "_o");
		}
	}
	
	n->setEvaluationFunction(&Workbench::sscore, this);
	n->setVerbose(false);
	n->setSilent(true);

	n->refine();
	
	double target = simpleScore(0.0);
	std::cout << "Parameters: " << n->parameterCount() <<
	"\tData/model: " << target << std::endl;
}

double Workbench::aveModel(Strain *strain, Serum *serum)
{
	Strain *orig = serum->strain();
	double strength = serum->strength();
	strength *= orig->averageEase();

	double offset = serum->offset() + orig->offset();
	double result = orig->aveVectorCompare(strain);

	double val = -result * strength;
	val += offset;

	return val;
}

double Workbench::modelForPair(Strain *strain, Serum *serum)
{
	Strain *orig = serum->strain();
	double strength = serum->strength();
	strength += orig->ease() - strain->ease();

	double offset = serum->offset() + orig->offset();
	double result = orig->vectorCompare(strain);
	result *= result;
	double val = -result * strength;

	val += offset;

	return val;
}

double Workbench::vectorModelForPair(Strain *strain, Serum *serum)
{
	return 0;
}

double Workbench::simpleScore(double weight)
{
	double diffs = 0;
	double count = 0;
	
	for (size_t i = 0; i < _tables.size(); i++)
	{
		for (size_t k = 0; k < _tables[i]->challengeCount(); k++)
		{
			Strain *challenge;
			Serum *serum;

			bool two = false;
			double data = _tables[i]->challenge(k, &challenge, &serum, &two);
			if (data != data)
			{
				continue;
			}
			
			if (challenge->isFree())
			{
				continue;
			}

			double val = modelForPair(challenge, serum);

			double scale = _tables[i]->scale();
			val *= scale;

			double diff = (val - data) * (val - data);
			
			if (!two && val > data)
			{
				diff = 0;
			}

			diffs += diff;

			count++;
		}
	}
	
	diffs /= count;

	return diffs;
}

double Workbench::gradientRefresh(void *object)
{
	Workbench *l = static_cast<Workbench *>(object);

	for (size_t i = 0; i < l->_strains.size(); i++)
	{
		l->_strains[i]->findPosition();
	}

	return 0;
}

double Workbench::vectorScore()
{
	return 0;
}

double Workbench::score()
{
	return simpleScore(0.00);
}

void Workbench::populateNames(char ****ptr)
{
	*ptr = (char ***)malloc(_sera.size() * sizeof(char ***));

	for (size_t i = 0; i < _sera.size(); i++)
	{
		(*ptr)[i] = (char **)malloc(_strains.size() * sizeof(char **));

		for (size_t j = 0; j < _strains.size(); j++)
		{
			std::string name = _sera[i]->name() + " / " + _strains[j]->name();

			(*ptr)[i][j] = (char *)malloc((name.length() + 1) * sizeof(char *));
			memcpy((*ptr)[i][j], &name[0], name.length() + 1);
		
		}
	}
}

void Workbench::populateRaw(double ***ptr, int type)
{
	*ptr = (double **)malloc(_sera.size() * sizeof(double **));

	for (size_t i = 0; i < _sera.size(); i++)
	{
		(*ptr)[i] = (double *)malloc(_strains.size() * sizeof(double));
		
		for (size_t j = 0; j < _strains.size(); j++)
		{
			if (!_strains[j]->hasSerum(_sera[i]))
			{
				(*ptr)[i][j] = NAN;
				continue;
			}

			double remove = _sera[i]->strain()->defaultOffset();
			
			double data = _strains[j]->serumValue(_sera[i]);
			data -= remove;

			if (type == 0)
			{
				(*ptr)[i][j] = data;
			}
			
			double scale = _strains[j]->serumScale(_sera[i]);
			double model = modelForPair(_strains[j], _sera[i]);
			model *= scale;
			model -= remove;
			
			if (type == 1)
			{
				(*ptr)[i][j] = model;
			}
			
			if (type == -1)
			{
				(*ptr)[i][j] = data - model;
			}

			(*ptr)[i][j] *= _scale;
			(*ptr)[i][j] += 0.5;
		}
	}
}

void Workbench::writeOut(std::string filename, int type)
{
	std::ofstream file;
	file.open(filename);
	file << ",";
	for (size_t j = 0; j < _strains.size(); j++)
	{
		file << _strains[j]->name() << ",";
	}
	file << std::endl;

	for (size_t i = 0; i < _sera.size(); i++)
	{
		file << _sera[i]->name() << ",";
		for (size_t j = 0; j < _strains.size(); j++)
		{
			if (!_strains[j]->hasSerum(_sera[i]))
			{
				file << "NAN,";
				continue;
			}

			double remove = _sera[i]->strain()->defaultOffset();
			double data = _strains[j]->serumValue(_sera[i]);
			data -= remove;

			double val = 0;

			if (type == 0)
			{
				val = data;
			}
			
			double model = modelForPair(_strains[j], _sera[i]);
			double scale = _strains[j]->serumScale(_sera[i]);
			model *= scale;
			model -= remove;
			
			if (type == 1)
			{
				val = model;
			}
			
			if (type == -1)
			{
				val = data - model;
			}

			val /= 2;

			file << val << ",";
		}
		file << std::endl;
	}
	
	file.close();
}

double Workbench::resultForVector(double *dir1, double *dir2)
{
	if (_scratch == NULL)
	{
		_scratch = new double[_dim];
	}

	for (size_t j = 0; j < _dim; j++)
	{
		double val = (dir1[j] - dir2[j]);
		_scratch[j] = val;
	}

	return resultForDirection(&_scratch[0]);
}

double Workbench::resultForDirection(double *dir)
{
	double pos = 0;
	double neg = 0;
	double sum = 0;

	for (size_t j = 0; j < _dim; j++)
	{
		double contrib = (dir[j] * dir[j]);
		sum += contrib;
		if (j >= _dim - 1)
		{
			neg += contrib;
		}
		else
		{
			pos += contrib;
		}
	}

	double sq = sqrt(sum);
	return sq;
}

double Workbench::resultForResidue(int i)
{
	return resultForDirection(_muts[i]->scalarPtr(0));
}

void Workbench::perResidueAntigenicity(std::string filename)
{
	std::ofstream file;
	file.open(filename);
	file << "mutation, ease, x, y, z, dist" << std::endl;
	for (size_t i = 0; i < _muts.size(); i++)
	{
		if (_muts[i]->silenced())
		{
			continue;
		}

		double sum = 0;
		file << _muts[i]->str() << ", ";
		file << _muts[i]->ease() << ", ";

		for (size_t j = 0; j < _dim; j++)
		{
//			double val = _muts[i]->scalar(j);
//			sum += val * val;
//			file << val << ", ";
		}
//		sum = sqrt(sum);
//		file << sum << ", " << std::endl;
	}
	file.close();
	
	file.open("distance_" + filename);

	for (size_t i = 0; i < _muts.size(); i++)
	{
		for (size_t j = 0; j < _muts.size(); j++)
		{
			double result = resultForVector(_muts[i]->scalarPtr(0),
			                                _muts[j]->scalarPtr(0));

			file << _muts[i]->str() << "," << _muts[j]->str() << "," << result << std::endl;
		}
	}

	file.close();
	
	file.open("strain_" + filename);
	file << "strain, ease, x, y, z" << std::endl;

	for (size_t i = 0; i < _strains.size(); i++)
	{
		_strains[i]->findPosition();
		vec3 dir = _strains[i]->aveDirection();
		
		file << _strains[i]->name() << ", ";
		file << _strains[i]->averageEase() << ", ";
		file << dir.x << ", " << dir.y << ", " << dir.z << ", ";
		file << std::endl;
	}

	file.close();
	
	file.open("strain_distance_" + filename);

	for (size_t i = 0; i < _strains.size() - 1; i++)
	{
		_strains[i]->findPosition();
		for (size_t j = i + 1; j < _strains.size(); j++)
		{
			_strains[j]->findPosition();
			double result = _strains[i]->vectorCompare(_strains[j]);
			if (result < 0) result = 0;

			file << _strains[i]->name() << "," << _strains[j]->name() << ",";
			file << result << std::endl;
		}
	}

	file.close();
}

void Workbench::antigenicity(std::string filename)
{
	std::ofstream file;
	file.open(filename);
	file << "mutation,ease,+-,antigenicity,+-" << std::endl;

	for (size_t i = 0; i < _muts.size(); i++)
	{
		double sum = _sums[i];
		double mean = sum / (double)_count;
		double sumsq = _sumsqs[i] / (double)_count;
		double stdev = sqrt(sumsq - mean * mean);
		
		if (_count > 1)
		{
			stdev /= sqrt(_count / (_count - 1));
		}
		else
		{
			stdev = 0;
		}
		
		double easedev = 0;
		double ease = _muts[i]->averageEase(&easedev);
		
		if (!_muts[i]->isUsed())
		{
			stdev = 0;
			mean = 0;
			easedev = 0;
			ease = 0;
		}

		file << _muts[i]->str() << ",";
		file << ease << "," << easedev << ",";
		file << mean << "," << stdev << std::endl;
	}

	{
		std::ofstream file;
		file.open("serum_strength.csv");
		file << "serum,strength" << std::endl;

		for (size_t i = 0; i < _sera.size(); i++)
		{
			file << _sera[i]->name() << ",";
			file << -_sera[i]->strain()->offset() + 
			_sera[i]->offset() << std::endl;
		}

		file.close();
	}
}

void Workbench::updateSums()
{
	_count++;
	
	double s = score();
	_scores.push_back(s);

	bool improved = false;
	if (s < _best)
	{
		_best = s;
		improved = true;
	}

	if (_sums.size() == 0)
	{
		_sums = std::vector<double>(_muts.size(), 0);
		_sumsqs = std::vector<double>(_muts.size(), 0);
	}

	for (size_t i = 0; i < _muts.size(); i++)
	{
		double result = resultForResidue(i);

		_sums[i] += result;
		_sumsqs[i] += result * result;
		_muts[i]->saveVector(-1, improved);
	}

	for (size_t i = 0; i < _strains.size(); i++)
	{
		_strains[i]->saveEase();
	}
	
	double ave = mean(_scores);
	double stdev = standard_deviation(_scores);
	std::cout << "Current scores: " << std::endl;
	std::cout << "\tTotal: " << _scores.size() << std::endl;
	std::cout << "\t Mean: " << ave << std::endl;
	std::cout << "\t Stdev: " << stdev << std::endl;
	
	
	rotations();

}

void Workbench::rotations()
{
	for (size_t i = 0; i < _muts[0]->trials(); i++)
	{
		for (size_t n = 0; n < _muts.size(); n++)
		{
			_muts[n]->loadSaved(i);
		}
		
		bestRotation();

		for (size_t n = 0; n < _muts.size(); n++)
		{
			_muts[n]->saveVector(i);
		}
	}
}

void Workbench::bestRotation()
{
	double *vals = (double *)calloc(_dim * _dim, sizeof(double));
	double **ptrs = (double **)malloc(sizeof(double *) * _dim);
	double *w = (double *)malloc(sizeof(double) * _dim);

	double *vVals = (double *)calloc(_dim * _dim, sizeof(double));
	double **vPtrs = (double **)malloc(sizeof(double *) * _dim);
	
	for (size_t i = 0; i < _dim; i++)
	{
		ptrs[i] = &vals[i * _dim];
		vPtrs[i] = &vVals[i * _dim];
	}

	for (size_t n = 0; n < _muts.size(); n++)
	{
		for (size_t j = 0; j < _dim; j++)
		{
			for (size_t i = 0; i < _dim; i++)
			{
				double add = _muts[n]->crossCorrelation(i, j);
				ptrs[i][j] += add;
			}
		}
	}

	int success = svdcmp((mat)ptrs, _dim, _dim, (vect) w, (mat) vPtrs);
	
	free(w);

	if (!success)
	{
		return;
	}

	double *rotvals = (double *)calloc(_dim * _dim, sizeof(double));
	double **rot = (double **)malloc(sizeof(double *) * _dim);
	
	for (size_t i = 0; i < _dim; i++)
	{
		rot[i] = &rotvals[i * _dim];
	}
	
	for (size_t j = 0; j < _dim; j++)
	{
		for (size_t i = 0; i < _dim; i++)
		{
			for (size_t n = 0; n < _dim; n++)
			{
				rot[i][j] += vPtrs[i][n] * ptrs[j][n];
			}
		}
	}
	
	for (size_t i = 0; i < _muts.size(); i++)
	{
		_muts[i]->rotateCurrent(rot);
	}
}

void Workbench::setup()
{
	_shouldSetup = false;
	reorderMutations();
	
	prepareVectors();
	degenerateSummary();
	startOffsets();
	prepare();
	
	gradientRefresh(this);
	score();

}

void Workbench::addMutations(std::vector<Mutation *> mutations)
{
	for (size_t i = 0; i < mutations.size(); i++)
	{
		if (std::find(_muts.begin(), _muts.end(),
		              mutations[i]) != _muts.end())
		{
			continue;
		}
		
		_muts.push_back(mutations[i]);
	}
	
	_shouldSetup = true;
}

void Workbench::addSera(std::vector<Serum *> sera)
{
	for (size_t i = 0; i < sera.size(); i++)
	{
		if (std::find(_sera.begin(), _sera.end(), sera[i]) != _sera.end())
		{
			continue;
		}
		
		_sera.push_back(sera[i]);
		_name2Serum[sera[i]->name()] = sera[i];
	}

	_shouldSetup = true;
}

void Workbench::addStrains(std::vector<Strain *> strains)
{
	for (size_t i = 0; i < strains.size(); i++)
	{
		if (std::find(_strains.begin(), _strains.end(), 
		              strains[i]) != _strains.end())
		{
			continue;
		}
		
		_strains.push_back(strains[i]);
		_name2Strain[strains[i]->name()] = strains[i];
	}

	_shouldSetup = true;
}

void Workbench::addTables(std::vector<Table *> tables)
{
	for (size_t i = 0; i < tables.size(); i++)
	{
		if (std::find(_tables.begin(), _tables.end(), 
		              tables[i]) != _tables.end())
		{
			continue;
		}
		
		_tables.push_back(tables[i]);
	}

	_shouldSetup = true;
}

Strain *Workbench::nextStrain(double *score)
{
	std::vector<Mutation *> copy = _muts;
	std::random_shuffle(copy.begin(), copy.end());
	Strain *tmp = new Strain("test");
	std::string list;
	for (size_t i = 0; i < 6; i++)
	{
		list += copy[i]->str() + ",";
	}
	list.pop_back();
	tmp->setList(list, _muts);

	double cc_sum = 0;
	double cc_with = 0;

	for (size_t i = 0; i < _muts.size() - 1; i++)
	{
		for (size_t j = i + 1; j < _muts.size(); j++)
		{
			double dots = 0;
			double amp_a = 0;
			double amp_b = 0;
			for (size_t k = 0; k < _strains.size(); k++)
			{
				double a = _strains[k]->hasMutation(_muts[i]);
				double b = _strains[k]->hasMutation(_muts[j]);

				dots += a * b;
				amp_a += a * a;
				amp_b += b * b;
			}

			double cc = dots / sqrt(amp_a * amp_b);
			cc_sum += cc;
			double a = tmp->hasMutation(_muts[i]);
			double b = tmp->hasMutation(_muts[j]);

			dots += a * b;
			amp_a += a * a;
			amp_b += a * a;

			double wcc = dots / sqrt(amp_a * amp_b);
			cc_with += wcc;
		}
	}
	
	*score = cc_with - cc_sum;
	return tmp;
}

void Workbench::clusterMutations(std::string filename)
{
	std::ofstream file;
	file.open(filename);
	
	std::map<Mutation *, std::map<Mutation*, double> > scores;
	
	for (size_t i = 0; i < _muts.size() - 1; i++)
	{
		for (size_t j = i + 1; j < _muts.size(); j++)
		{
			double dots = 0;
			double amp_a = 0;
			double amp_b = 0;
			for (size_t k = 0; k < _strains.size(); k++)
			{
				double a = _strains[k]->hasMutation(_muts[i]);
				double b = _strains[k]->hasMutation(_muts[j]);
				
				dots += a * b;
				amp_a += a * a;
				amp_b += b * b;
			}
			
			double cc = dots / sqrt(amp_a * amp_b);
			scores[_muts[i]][_muts[j]] = cc;
			scores[_muts[j]][_muts[i]] = cc;
		}
	}

	for (size_t i = 0; i < _muts.size() - 1; i++)
	{
		for (size_t j = i + 1; j < _muts.size(); j++)
		{
			double dots = 0;
			double amp_a = 0;
			double amp_b = 0;

			for (size_t k = 0; k < _muts.size(); k++)
			{
				double a = scores[_muts[i]][_muts[k]];
				double b = scores[_muts[j]][_muts[k]];
				dots += a * b;
				amp_a += a * a;
				amp_b += b * b;
			}
			double cc = dots / sqrt(amp_a * amp_b);

			file << _muts[i]->str() << "," << _muts[j]->str() << ","
			<< cc << std::endl;
		}
	}

	file.close();
	
	for (size_t j = 0; j < 5; j++)
	{
		double best = 0;
		Strain *best_strain = NULL;
		for (size_t i = 0; i < 10000; i++)
		{
			double score = 0;
			Strain *tmp = nextStrain(&score);
			if (score > best)
			{
				best_strain = tmp;
				best = score;
			}
		}
		
		std::cout << best_strain->summary() << "\t" << best << std::endl;
		_strains.push_back(best_strain);
	}

}

void Workbench::makeStrainFree(std::string str)
{
	to_lower(str);

	if (_name2Strain.count(str) == 0)
	{
		return;
	}

	Strain *strain = _name2Strain[str];
	strain->setFree(true);
	std::cout << "Set " << strain->name() << " to free set." << std::endl;
}

void Workbench::markStrains(std::string strains)
{

}

Settings Workbench::displaySettings(std::string filename)
{
	Settings settings(this, filename);
	if (settings.isValid())
	{
		settings.apply();
	}

	return settings;
}

std::vector<Strain *> Workbench::getStrains(std::string list)
{
	std::vector<std::string> words = split(list, ',');
	std::vector<Strain *> strains;

	for (size_t i = 0; i < words.size(); i++)
	{
		Strain *strain = _name2Strain[words[i]];
		if (strain)
		{
			strains.push_back(strain);
		}
	}
	
	return strains;
}


void Workbench::assignToPDB(std::string pdb)
{

}

