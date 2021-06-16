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

#include "Loader.h"
#include "Strain.h"
#include "Serum.h"
#include <fstream>
#include <algorithm>
#include <hcsrc/FileReader.h>
#include <hcsrc/RefinementNelderMead.h>
#include <hcsrc/RefinementLBFGS.h>
#include <hcsrc/Any.h>
#include <libsrc/shared_ptrs.h>

int Loader::_dim = 2;

Loader::Loader()
{
	_type = ModelVector;
	_refineOffset = false;
	_refineStrength = false;
	_count = 0;
	_scale = 1;
}

void Loader::load(std::string filename)
{
	if (!file_exists(filename))
	{
		std::cout << "Could not find " << filename << std::endl;
		return;
	}

	std::string contents = get_file_contents(filename);
	std::vector<std::string> lines = split(contents, '\n');
	std::cout << std::endl;

	for (size_t i = 0; i < lines.size(); i++)
	{
		std::string line = lines[i];
		to_lower(line);

		if (line.rfind("strain ", 0) == 0)
		{
			line.erase(0, strlen("strain "));
			defineStrain(line);
		}

		if (line.rfind("serum ", 0) == 0)
		{
			line.erase(0, strlen("serum "));
			defineSerum(line);
		}

		if (line.rfind("challenge ", 0) == 0)
		{
			line.erase(0, strlen("challenge "));
			defineChallenge(line);
		}

		if (line.rfind("cooperative ", 0) == 0)
		{
			line.erase(0, strlen("cooperative "));
			enablePairs(line);
		}
	}
	
	reorderMutations();
	std::cout << std::endl;
	
	fillInZeros();
	prepareVectors();
	degenerateSummary();
	prepare();
	score();
}

void Loader::fillInZeros()
{
	for (size_t i = 0; i < _sera.size(); i++)
	{
		Strain *strain = _sera[i]->strain();
		_strainMap[strain][_sera[i]] = 0;
		_freeMap[strain][_sera[i]] = 0;
	}
}

void Loader::removeUnused()
{
	for (size_t i = 0; i < _muts.size(); i++)
	{
		bool found = false;

		if (_muts[i].find('+') == std::string::npos)
		{
			continue;
		}
		
		for (size_t k = 0; k < _sera.size() && !found; k++)
		{
			for (size_t j = 0; j < _strains.size(); j++)
			{
				std::vector<std::string> all;
				Strain *a = _sera[k]->strain();
				Strain *b = _strains[j];

				if (Strain::hasCooperative(a, b, _muts[i]))
				{
					found = true;
					break;
				}
			}
		}

		if (!found)
		{
			_muts.erase(_muts.begin() + i);
			i--;
		}
	}
	
	prepareVectors();
	degenerateSummary();
}

void Loader::prepareVectors()
{
	std::cout << "Vectors: " << std::endl;
	
	for (size_t i = 0; i < _strains.size(); i++)
	{
		_strains[i]->clearPrecalculated();
		std::vector<double> vecs = _strains[i]->generateVector(_muts);
	}

	std::cout << std::endl;
	
	std::cout << "All " << _muts.size() << " mutations over " << _strains.size() 
	<< " strains: " << std::endl;

	for (size_t i = 0; i < _muts.size(); i++)
	{
		std::cout << _muts[i] << " " << std::flush;
	}

	std::cout << std::endl;

	size_t count = 0;
	
	for (size_t i = 0; i < _strains.size(); i++)
	{
		for (size_t j = 0; j < _sera.size(); j++)
		{
			if (_strainMap.count(_strains[i]) &&
			    _strainMap[_strains[i]].count(_sera[j]))
			{
				count++;
			}
		}
	}
	std::cout << "Total number of observations: " << count << std::endl;
}

void Loader::degenerateSummary()
{
	std::map<std::string, std::vector<std::string>> map;
	for (size_t i = 0; i < _muts.size(); i++)
	{
		std::string m = _muts[i];
		std::string v;
		
		for (size_t j = 0; j < _strains.size(); j++)
		{
			bool has = (_strains[j]->hasMutation(m));
			v += (has ? "1" : "0");
		}

		map[v].push_back(m);
	}
	
	std::map<std::string, std::vector<std::string>>::iterator it;

	std::cout << std::endl;
	std::cout << "Mutations grouped by degeneracy:" << std::endl;
	for (it = map.begin(); it != map.end(); it++)
	{
		std::vector<std::string> v = it->second;
		
		for (size_t j = 0; j < v.size(); j++)
		{
			std::cout << v[j] << ", ";
		}

		std::cout << std::endl;
	}
	
	std::cout << std::endl;
}

typedef struct
{
	std::string str;
	int res;
} StrRes;

bool higher_seq(StrRes &a, StrRes &b)
{
	if (a.str.find('+') != std::string::npos && 
	    b.str.find('+') == std::string::npos)
	{
		return false;
	}
	if (b.str.find('+') != std::string::npos && 
	    a.str.find('+') == std::string::npos)
	{
		return true;
	}
	return (a.res < b.res);
}

void Loader::reorderMutations()
{
	std::vector<StrRes> residues;
	for (size_t i = 0; i < _muts.size(); i++)
	{
		StrRes res;
		res.str = _muts[i];
		char *old = &_muts[i][0];
		char *pos = old;
		int val = 0;
		
		while (pos == old)
		{
			val = strtol(old, &pos, 10);
			pos++;
			old++;
		}
		
		res.res = val;
		residues.push_back(res);
	}

	std::sort(residues.begin(), residues.end(), higher_seq);
	
	for (size_t i = 0; i < _muts.size(); i++)
	{
		_muts[i] = residues[i].str;
	}
}

void Loader::defineSerum(std::string str)
{
	trim(str);
	std::vector<std::string> equation = split(str, '=');
	if (equation.size() != 2)
	{
		std::cout << "Strange serum definition: " << str << std::endl;
		return;
	}
	
	trim(equation[0]);
	trim(equation[1]);
	Strain *one = _name2Strain[equation[1]];
	if (one == NULL)
	{
		std::cout << "WARNING: Could not find " << equation[1] 
		<< ", ignoring." << std::endl;
		return;
	}
	
	Serum *serum = new Serum(equation[0], one);
	_sera.push_back(serum);
	_name2Serum[serum->name()] = serum;

	std::cout << "Defining " << serum->name() << " as " <<
	one->name() << std::endl;
}

void Loader::defineStrain(std::string str)
{
	trim(str);
	std::vector<std::string> equation = split(str, '=');
	if (str.back() == '=')
	{
		equation.push_back("");
		
	}
	if (equation.size() != 2)
	{
		std::cout << "Strange strain equation: " << str << std::endl;
		return;
	}
	
	trim(equation[0]);
	trim(equation[1]);
	
	for (size_t i = 0; i < _strains.size(); i++)
	{
		if (_strains[i]->name() == equation[0])
		{
			return;
		}
	}

	Strain *strain = new Strain(equation[0]);
	strain->setList(equation[1]);

	_strains.push_back(strain);
	_name2Strain[strain->name()] = strain;

	strain->addToMuts(_muts);
}

void Loader::defineChallenge(std::string str)
{
	std::vector<std::string> challenges = split(str, ',');
	if (challenges.size() < 3)
	{
		std::cout << "Strange challenge: " << str << std::endl;
		return;
	}
	
	trim(challenges[0]);
	trim(challenges[1]);
	trim(challenges[2]);

	Strain *one = _name2Strain[challenges[0]];
	Serum *two = _name2Serum[challenges[1]];
	if (one == NULL)
	{
		std::cout << "WARNING: Could not find " << challenges[0] 
		<< ", ignoring." << std::endl;
		return;
	}
	if (two == NULL)
	{
		std::cout << "WARNING: Could not find " << challenges[1]
		<< ", ignoring." << std::endl;
		return;
	}

	double val = atof(challenges[2].c_str());
	
	val /= _scale;

	_strainMap[one][two] = val;
	_freeMap[one][two] = false;
	
	if (challenges.size() >= 4)
	{
		std::string ch = challenges[3];
		trim(ch);
		if (ch == "free")
		{
			_freeMap[one][two] = true;
//			std::cout << "(free) ";
		}
	}

//	std::cout << one->name() << " strain [v] " << two->name() << 
//	" serum -> " << val << std::endl;
}

void Loader::refineLoop(int count)
{
	std::cout << "*** REFINE ***" << std::endl;
	srand(time(NULL));
	_reals.clear();

	for (size_t i = 0; i < count; i++)
	{
		refine();
	}
	
	updateSums();
	
	addResult();

	emit resultReady();
}

void Loader::addResult()
{
	std::string add = "";
	add += "result_" + f_to_str(score(), 3) + "_";
	add += i_to_str(_count) + ",";
	
	for (size_t i = 0; i < _reals.size(); i++)
	{
		double result = resultForResidue(i);
		add += f_to_str(result, 3) + ",";
	}

	add += "\n";

	_results += add;
}

void Loader::writeResultVectors(std::string filename)
{
	std::ofstream file;
	file.open(filename);
	
	file << "result_num,";
	for (size_t i = 0; i < _muts.size(); i++)
	{
		file << _muts[i] << ",";
	}
	file << std::endl;
	
	file << _results;
	file.close();
}

void Loader::prepare()
{
	if (_reals.size() == 0 || _reals.size() != _muts.size())
	{
		_reals = std::vector<double>(_muts.size(), 0.2);
		_scales = std::vector<std::vector<double> >(_muts.size(), 
		                                            std::vector<double>(_dim, 0));

		for (size_t i = 0; i < _reals.size(); i++)
		{
			_reals[i] = ((rand() / (double)RAND_MAX) - 0.5) / 4;

			for (size_t j = 0; j < _dim; j++)
			{
				_scales[i][j] = ((rand() / (double)RAND_MAX) - 0.5) / 4;
				if (j >= 2)
				{
					_scales[i][j] /= 10;
				}
			}
		}
		
		for (size_t i = 0; i < _sera.size(); i++)
		{
			(*_sera[i]->strengthPtr()) = 1;
			(*_sera[i]->offsetPtr()) = 0;
		}
	}
}

void Loader::refine()
{
	prepare();

	double s = score();
	RefinementLBFGSPtr n = RefinementLBFGSPtr(new RefinementLBFGS());
	n->setGradientRefresh(this, gradientRefresh);
	n->setCycles(1000);
	
	for (size_t i = 0; i < _anys.size(); i++)
	{
		delete _anys[i];
	}
	_anys.clear();

	for (size_t i = 0; i < _reals.size(); i++)
	{
		if (_type == ModelSimple)
		{
			Any *any_real = new Any(&_reals[i]);
			_anys.push_back(any_real);

			std::string mut = _muts[i];

			n->addParameter(any_real, Any::get, Any::set, 0.2, 0.02, mut + "_r");
		}

		if (_type == ModelVector)
		{
			for (size_t j = 0; j < _dim; j++)
			{
				Any *any_real = new Any(&_scales[i][j]);
				_anys.push_back(any_real);

				std::string mut = _muts[i];

				n->addParameter(any_real, Any::get, Any::set, 0.2, 0.02, 
				                mut + "_" + i_to_str(j));
			}
		}

	}
	
	for (size_t i = 0; i < _sera.size(); i++)
	{
		std::string name = _sera[i]->name();

		if (_refineStrength)
		{
			double *strength = _sera[i]->strengthPtr();
			Any *any = new Any(strength);
			_anys.push_back(any);

			n->addParameter(any, Any::get, Any::set, 0.5, 0.05, name + "_s");
		}

		if (_refineOffset)
		{
			double *offset = _sera[i]->offsetPtr();
			Any *any = new Any(offset);
			_anys.push_back(any);

			n->addParameter(any, Any::get, Any::set, 0.2, 0.02, name + "_o");
		}
	}
	
	n->setEvaluationFunction(&Loader::sscore, this);
	n->setVerbose(false);
	n->setSilent(true);

	n->refine();
	
	double target = simpleScore(0.0);
	double both = score();
	std::cout << "Parameters: " << n->parameterCount() <<
	"\tData/model: " << target
	<< "\trestraints alone: " << both - target
	<< "\tboth: " << both << std::endl;
}

double Loader::modelForPair(Strain *strain, Serum *serum)
{
	Strain *orig = serum->strain();
	double strength = serum->strength() / 1;
	double offset = serum->offset();
	double result = 0;
	
	if (_type == ModelSimple)
	{
		result = strain->compareToStrain(orig, _reals, _muts);
	}
	else
	{
		result = strain->vectorCompare(orig, _scales);
	}

	double val = -result * strength;
	val -= offset;

	return val;
}

double Loader::vectorModelForPair(Strain *strain, Serum *serum)
{
	return 0;
}

double Loader::simpleScore(double weight)
{
	double diffs = 0;
	double count = 0;
	double restraints = 0;
	
	gradientRefresh(this);
	
	for (size_t i = 0; i < _strains.size(); i++)
	{
		Strain *challenge = _strains[i];
		for (size_t j = 0; j < _sera.size(); j++)
		{
			Serum *serum = _sera[j];
			Strain *orig = serum->strain();
			
			if (challenge == orig)
			{
				continue;
			}
			
			if (!_strainMap.count(challenge) ||
			    !_strainMap[challenge].count(serum)
			    || _freeMap[challenge][serum])
			{
				continue;
			}

			double data = _strainMap[challenge][serum];
			double val = modelForPair(challenge, serum);

			double diff = (val - data) * (val - data);
			diffs += diff;
			
			count++;
		}
	}
	
	for (size_t i = 0; i < _muts.size(); i++)
	{
		restraints += restraintForResidue(i);
	}

	diffs /= count;
	restraints /= count;
	restraints = sqrt(restraints);

	return diffs + restraints * weight;
}

double Loader::gradientRefresh(void *object)
{
	Loader *l = static_cast<Loader *>(object);

	for (size_t i = 0; i < l->_strains.size() && l->_type == ModelVector; i++)
	{
		l->_strains[i]->findPosition(l->_scales);
	}

	return 0;
}

double Loader::vectorScore()
{
	return 0;
}

double Loader::score()
{
	return simpleScore(0.00);
}

void Loader::populateNames(char ****ptr)
{
	*ptr = (char ***)malloc(_sera.size() * sizeof(char ***));

	for (size_t i = 0; i < _sera.size(); i++)
	{
		(*ptr)[i] = (char **)malloc(_strains.size() * sizeof(char **));

		for (size_t j = 0; j < _strains.size(); j++)
		{
			std::string name = _sera[i]->name() + " / " + _strains[j]->name();

			(*ptr)[i][j] = (char *)malloc(name.length() * sizeof(char *));
			memcpy((*ptr)[i][j], &name[0], name.length());
		
		}
	}
}

void Loader::populateRaw(double ***ptr, int type)
{
	*ptr = (double **)malloc(_sera.size() * sizeof(double **));

	for (size_t i = 0; i < _sera.size(); i++)
	{
		(*ptr)[i] = (double *)malloc(_strains.size() * sizeof(double));
		
		for (size_t j = 0; j < _strains.size(); j++)
		{
			if (!_strainMap.count(_strains[j]) ||
			    !_strainMap[_strains[j]].count(_sera[i]))
			{
				(*ptr)[i][j] = NAN;
				continue;
			}

			double data = _strainMap[_strains[j]][_sera[i]];

			if (type == 0)
			{
				(*ptr)[i][j] = data;
			}
			
			double model = modelForPair(_strains[j], _sera[i]);
			
			if (type == 1)
			{
				(*ptr)[i][j] = model;
			}
			
			if (type == -1)
			{
				(*ptr)[i][j] = data - model;
			}

			(*ptr)[i][j] /= 2;
			(*ptr)[i][j] += 0.5;
		}
	}
}

void Loader::writeOut(std::string filename, int type)
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
			if (!_strainMap.count(_strains[j]) ||
			    !_strainMap[_strains[j]].count(_sera[i]))
			{
				file << "NAN,";
				continue;
			}

			double data = _strainMap[_strains[j]][_sera[i]];
			double val = 0;

			if (type == 0)
			{
				val = data;
			}
			
			double model = modelForPair(_strains[j], _sera[i]);
			
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

void Loader::setModelToAverage()
{
	for (size_t i = 0; i < _sums.size(); i++)
	{
		double sum = _sums[i];
		double mean = sum / (double)_count;
		_reals[i] = mean;
	}
}

double Loader::restraintForResidue(int i)
{
	double sum = 0;
	for (size_t j = 0; j < _dim; j++)
	{
		double contrib = (_scales[i][j] * _scales[i][j]);
		sum += contrib;
	}
	
	return sum * sum;
}

double Loader::resultForVector(std::vector<double> dir1, 
                               std::vector<double> dir2)
{
	std::vector<double> subtract;

	for (size_t j = 0; j < _dim; j++)
	{
		subtract.push_back(dir1[j] - dir2[j]);
	}

	return resultForDirection(subtract);
}

double Loader::resultForDirection(std::vector<double> dir)
{
	double sum = 0;
	double pos = 0;
	double neg = 0;

	for (size_t j = 0; j < _dim; j++)
	{
		double contrib = (dir[j] * dir[j]);
		sum += contrib;
		if (j >= 2)
		{
			neg += contrib;
		}
		else
		{
			pos += contrib;
		}
	}

	return (pos - neg);
}

double Loader::resultForResidue(int i)
{
	if (_type == ModelSimple)
	{
		return _reals[i];
	}
	
	return resultForDirection(_scales[i]);
}

void Loader::perResidueAntigenicity(std::string filename)
{
	if (_type == ModelSimple)
	{
		return;
	}

	std::ofstream file;
	file.open(filename);
	file << "mutation, x, y, z" << std::endl;
	for (size_t i = 0; i < _scales.size(); i++)
	{
		file << _muts[i] << ", ";
		for (size_t j = 0; j < _dim; j++)
		{
			file << _scales[i][j] << ", ";
		}
		file << std::endl;
	}
	file.close();
	
	file.open("distance_" + filename);

	for (size_t i = 0; i < _muts.size(); i++)
	{
		for (size_t j = 0; j < _muts.size(); j++)
		{
			double result = resultForVector(_scales[i], _scales[j]);

			file << _muts[i] << "," << _muts[j] << "," << result << std::endl;
		}
	}

	file.close();
	
	file.open("strain_" + filename);
	file << "strain, x, y, z" << std::endl;

	for (size_t i = 0; i < _strains.size(); i++)
	{
		_strains[i]->findPosition(_scales);
		std::vector<double> dir = _strains[i]->direction();
		
		file << _strains[i]->name() << ", ";
		for (size_t j = 0; j < dir.size(); j++)
		{
			file << dir[j] << ", ";
		}
		file << std::endl;
	}

	file.close();
	
	file.open("strain_distance_" + filename);

	for (size_t i = 0; i < _strains.size(); i++)
	{
		_strains[i]->findPosition(_scales);
		for (size_t j = 0; j < _strains.size(); j++)
		{
			_strains[j]->findPosition(_scales);
			double result = _strains[i]->vectorCompare(_strains[j], _scales);
			if (result < 0) result = 0;

			file << _strains[i]->name() << "," << _strains[j]->name() << ",";
			file << result << std::endl;
		}
	}

	file.close();
}

void Loader::antigenicity(std::string filename)
{
	std::ofstream file;
	file.open(filename);
	file << "mutation,antigenicity,+-" << std::endl;

	for (size_t i = 0; i < _sums.size(); i++)
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

		file << _muts[i] << "," << mean << "," << stdev << std::endl;
	}
}

void Loader::updateSums()
{
	_count++;

	if (_sums.size() == 0)
	{
		_sums = std::vector<double>(_reals.size(), 0);
		_sumScales = std::vector<std::vector<double>>(_reals.size(), 
		                                              std::vector<double>(_dim, 
		                                              0));
		_sumScaleSqs = std::vector<std::vector<double>>(_reals.size(), 
		                                                std::vector<double>(_dim, 
		                                                0));
		_sumsqs = std::vector<double>(_reals.size(), 0);
	}

	for (size_t i = 0; i < _reals.size(); i++)
	{
		double result = resultForResidue(i);
		
		for (size_t j = 0; j < _dim; j++)
		{
			_sumScales[i][j] += _scales[i][j];
			_sumScaleSqs[i][j] += _scales[i][j] * _scales[i][j];
		}

		_sums[i] += result;
		_sumsqs[i] += result * result;
	}
}

void Loader::enablePairs(std::string line)
{
	std::vector<std::string> muts = split(line, ',');
	std::cout << "Enabling cooperation between mutations ";
	
	for (size_t i = 0; i < muts.size(); i++)
	{
		std::string tmp = muts[i];
		trim(tmp);
		muts[i] = tmp;
	}
	
	std::sort(muts.begin(), muts.end(), std::greater<std::string>());

	for (size_t i = 0; i < muts.size() - 1; i++)
	{
		std::cout << muts[i] << ", ";
		for (size_t j = i + 1; j < muts.size(); j++)
		{
			std::string combo = muts[i] + "+" + muts[j];
			_muts.push_back(combo);
		}
	}
	
	std::cout << "and " << muts.back() << "." <<  std::endl;
}
