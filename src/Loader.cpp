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
#include "Mutation.h"
#include "Serum.h"
#include <fstream>
#include <algorithm>
#include <hcsrc/FileReader.h>
#include <hcsrc/RefinementNelderMead.h>
#include <hcsrc/RefinementLBFGS.h>
#include <hcsrc/Any.h>
#include <libsrc/shared_ptrs.h>

int Loader::_dim = 2;
double *Loader::_scratch = NULL;

Loader::Loader()
{
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
	}
	
	reorderMutations();
	std::cout << std::endl;
	
	fillInZeros();
	prepareVectors();
	degenerateSummary();
	prepare();
	
	gradientRefresh(this);
	score();
}

void Loader::fillInZeros()
{
	for (size_t i = 0; i < _sera.size(); i++)
	{
		Strain *strain = _sera[i]->strain();
		
		if (strain->hasSerum(_sera[i]))
		{
			continue;
		}
		
		strain->addSerum(_sera[i], 0, false);
	}
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
		std::cout << _muts[i]->str() << " " << std::flush;
	}

	std::cout << std::endl;

	size_t count = 0;
	
	for (size_t i = 0; i < _strains.size(); i++)
	{
		count += serumCount();
	}

	std::cout << "Total number of observations: " << count << std::endl;
}

void Loader::degenerateSummary()
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
	return (a.res < b.res);
}

void Loader::reorderMutations()
{
	std::vector<StrRes> residues;
	for (size_t i = 0; i < _muts.size(); i++)
	{
		StrRes res;
		res.mut = _muts[i];
		std::string mut = _muts[i]->str();
		char *old = &mut[0];
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
		_muts[i] = residues[i].mut;
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
	strain->setLoader(this);
	strain->setList(equation[1], _muts);

	_strains.push_back(strain);
	_name2Strain[strain->name()] = strain;
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
	bool free = false;

	if (challenges.size() >= 4)
	{
		std::string ch = challenges[3];
		trim(ch);
		if (ch == "free")
		{
			free = true;
		}
	}

	one->addSerum(two, val, free);
}

void Loader::refineLoop(int count)
{
	std::cout << "*** REFINE ***" << std::endl;
	srand(time(NULL));

	prepare();

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
	
	for (size_t i = 0; i < _muts.size(); i++)
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
		file << _muts[i]->str() << ",";
	}
	file << std::endl;
	
	file << _results;
	file.close();
}

void Loader::prepare()
{
	for (size_t i = 0; i < _muts.size(); i++)
	{
		_muts[i]->randomiseVector();
	}

	for (size_t i = 0; i < _sera.size(); i++)
	{
		(*_sera[i]->strengthPtr()) = 1;
		(*_sera[i]->offsetPtr()) = 0;
	}

	std::cout << "Starting score: " << score() << std::endl;
}

void Loader::refine()
{
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
		for (size_t j = 0; j < _dim; j++)
		{
			Any *any_real = new Any(_muts[i]->scalarPtr(j));
			_anys.push_back(any_real);
			any_real->setRefresh(Mutation::refresh, _muts[i]);

			std::string mut = _muts[i]->str();

			n->addParameter(any_real, Any::get, Any::set, 0.2, 0.02, 
			                mut + "_" + i_to_str(j));
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
			any->setRefresh(Strain::staticFindPosition, _sera[i]->strain());

			n->addParameter(any, Any::get, Any::set, 0.5, 0.05, name + "_s");
		}

		if (_refineOffset)
		{
			double *offset = _sera[i]->offsetPtr();
			Any *any = new Any(offset);
			_anys.push_back(any);
			any->setRefresh(Strain::staticFindPosition, _sera[i]->strain());

			n->addParameter(any, Any::get, Any::set, 0.2, 0.02, name + "_o");
		}
	}
	
	n->setEvaluationFunction(&Loader::sscore, this);
	n->setVerbose(false);
	n->setSilent(true);

	n->refine();
	
	double target = simpleScore(0.0);
	std::cout << "Parameters: " << n->parameterCount() <<
	"\tData/model: " << target << std::endl;
}

double Loader::modelForPair(Strain *strain, Serum *serum)
{
	Strain *orig = serum->strain();
	double strength = serum->strength() / 1;
	double offset = serum->offset();
	double result = strain->vectorCompare(orig);

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
	
	gradientRefresh(this);
	
	for (size_t i = 0; i < _strains.size(); i++)
	{
		Strain *challenge = _strains[i];
		for (size_t j = 0; j < challenge->serumCount(); j++)
		{
			Serum *serum = challenge->serum(j);
			Strain *orig = serum->strain();
			
			if (challenge == orig)
			{
//				continue;
			}
			
			if (challenge->freeSerum(serum))
			{
				continue;
			}

			double data = challenge->serumValue(serum);
			double val = modelForPair(challenge, serum);

			double diff = (val - data) * (val - data);
			diffs += diff;
			
			count++;
		}
	}

	diffs /= count;

	return diffs;
}

double Loader::gradientRefresh(void *object)
{
	Loader *l = static_cast<Loader *>(object);

	for (size_t i = 0; i < l->_strains.size(); i++)
	{
		l->_strains[i]->findPosition();
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
			if (!_strains[j]->hasSerum(_sera[i]))
			{
				(*ptr)[i][j] = NAN;
				continue;
			}

			double data = _strains[j]->serumValue(_sera[i]);

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
			if (!_strains[j]->hasSerum(_sera[i]))
			{
				file << "NAN,";
				continue;
			}

			double data = _strains[j]->serumValue(_sera[i]);
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

double Loader::resultForVector(double *dir1, double *dir2)
{
	if (_scratch == NULL)
	{
		_scratch = new double[_dim];
	}

	for (size_t j = 0; j < _dim; j++)
	{
		_scratch[j] = dir1[j] - dir2[j];
	}

	return resultForDirection(&_scratch[0]);
}

double Loader::resultForDirection(double *dir)
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

	double sq = sqrt(pos) - sqrt(neg);
	return sq;
}

double Loader::resultForResidue(int i)
{
	return resultForDirection(_muts[i]->scalarPtr(0));
}

void Loader::perResidueAntigenicity(std::string filename)
{
	std::ofstream file;
	file.open(filename);
	file << "mutation, x, y, z" << std::endl;
	for (size_t i = 0; i < _muts.size(); i++)
	{
		file << _muts[i]->str() << ", ";
		for (size_t j = 0; j < _dim; j++)
		{
			file << _muts[i]->scalar(j) << ", ";
		}
		file << std::endl;
	}
	file.close();
	
	file.open("distance_" + filename);

	for (size_t i = 0; i < _muts.size(); i++)
	{
		for (size_t j = 0; j < _muts.size(); j++)
		{
			double result = resultForVector(_muts[i]->scalarPtr(i),
			                                _muts[j]->scalarPtr(j));

			file << _muts[i]->str() << "," << _muts[j]->str() << "," << result << std::endl;
		}
	}

	file.close();
	
	file.open("strain_" + filename);
	file << "strain, x, y, z" << std::endl;

	for (size_t i = 0; i < _strains.size(); i++)
	{
		_strains[i]->findPosition();
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
		_strains[i]->findPosition();
		for (size_t j = 0; j < _strains.size(); j++)
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

void Loader::antigenicity(std::string filename)
{
	std::ofstream file;
	file.open(filename);
	file << "mutation,antigenicity,+-" << std::endl;

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

		file << _muts[i]->str() << "," << mean << "," << stdev << std::endl;
	}
}

void Loader::updateSums()
{
	_count++;

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
	}
}
