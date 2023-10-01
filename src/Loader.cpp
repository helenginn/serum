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
#include "Table.h"
#include "Strain.h"
#include "Mutation.h"
#include "Workbench.h"
#include "Serum.h"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <hcsrc/FileReader.h>

Loader::Loader(Workbench *b)
{
	_bench = b;
	_forgiving = true;
	_currentTable = NULL;
}

void Loader::parseSpreadsheet(std::string folder)
{
	_folder = folder;
	_currentTable = new Table(folder);
	std::string search = folder + "/*.csv";
	std::vector<std::string> csvs;

	try
	{
		csvs = glob_pattern(search);
	}
	catch (std::runtime_error &err)
	{
		std::cout << err.what();
		return;
	}
	
	_forgiving = true;

	for (size_t i = 0; i < csvs.size(); i++)
	{
		std::cout << csvs[i] << std::endl;
		loadCSV(csvs[i]);
	}

	_forgiving = false;

	for (size_t i = 0; i < csvs.size(); i++)
	{
		loadCSV(csvs[i]);
	}
	
	_tables.push_back(_currentTable);
	_currentTable = NULL;
	
	std::cout << _sera.size() << " sera" << std::endl;
	std::cout << _strains.size() << " strains" << std::endl;
	std::cout << _muts.size() << " muts" << std::endl;

	_bench->addMutations(_muts);
	_bench->addSera(_sera);
	_bench->addStrains(_strains);
	_bench->addTables(_tables);
}

void Loader::loadCSV(std::string &filename)
{
	_filename = filename;
	std::string contents = get_file_contents(filename);
	std::string fixed = defenestrate(contents);
	to_lower(fixed);
	std::vector<std::string> lines = split(fixed, '\n');

	size_t result = 0;
	while (result < lines.size())
	{
		result = doLine(lines, result);
	}
}

size_t Loader::doChallenge(std::vector<std::string > &lines, size_t line)
{
	std::vector<std::string> bits = split(lines[line], ',');
	std::vector<Strain *> strainList;
	
	for (size_t i = 1; i < bits.size(); i++)
	{
		trim(bits[i]);
		
		if (bits[i].length() == 0 || bits[i][0] == '#')
		{
			strainList.push_back(NULL);
			continue;
		}

		if (_name2Strain.count(bits[i]) == 0)
		{
			std::string error = "Strain requested for challenge table," +
			bits[i] + ", does not exist!";
			throw error;
		}

		Strain *str = _name2Strain[bits[i]];
		
		strainList.push_back(str);
	}

	for (size_t i = 0; i < strainList.size(); i++)
	{
		if (strainList[i] == NULL)
		{
			continue;
		}

		std::cout << strainList[i]->name() << std::endl;
	}
	std::cout << std::endl;
	
	size_t next = line + 1;
	
	while (next < lines.size())
	{
		lines[next] += ',';
		bits = split(lines[next], ',');

		if (bits.size() == 0 || (bits[0].length() > 0 && bits[0][0] == '#'))
		{
			next++;
			continue;
		}
		
		if (bits[0] == "strain" || bits[1] == "challenge")
		{
			break;
		}

		bits.resize(strainList.size() + 1);
		trim(bits[0]);
		
		bool blank = true;
		for (size_t j = 0; j < bits.size(); j++)
		{
			trim(bits[j]);
			if (bits[j].length() > 0)
			{
				blank = false;
			}
		}
		
		if (blank)
		{
			next++;
			continue;
		}
		
		if (_name2Serum.count(bits[0]) == 0)
		{
			std::string error = "Serum requested for challenge table," +
			bits[0] + ", does not exist!";
			throw error;
		}

		Serum *serum = _name2Serum[bits[0]];
		
		for (size_t j = 1; j < bits.size(); j++)
		{
			trim(bits[j]);
			if (bits[j].length() == 0)
			{
				continue;
			}

			int idx = j - 1;
			Strain *chosen = strainList[idx];
			if (chosen == NULL)
			{
				/*
				std::string error = "Serum/strain value given for " +
				bits[0] + " has missing header!";
				std::cout << error << std::endl;
				throw -1;
				*/
				continue;
			}
			
			int dilution = atof(bits[j].c_str());
			double val = log(dilution) / log(10);
			bool two = (dilution > 1e-6);
			
			if (dilution <= 1)
			{
				val = log(1) / log(10);
			}

			_currentTable->addValue(chosen, serum, val, two);
		}
		
		next++;
	}

	return next + 1;
}

size_t Loader::doLine(std::vector<std::string> &lines, size_t line)
{
	std::vector<std::string> bits = split(lines[line], ',');
	
	if (bits.size() == 0 || (bits[0].length() > 0 && bits[0][0] == '#'))
	{
		return line + 1;
	}
	
	size_t next = line + 1;

	if (_forgiving && bits[0].rfind("strain", 0) == 0)
	{
		try
		{
			next = doStrain(lines, line);
		}
		catch (int e)
		{
			std::cout << "Unrecoverable error" << " ";
			std::cout << _filename << " : " << next << std::endl;
			exit(0);
		}
	}

	bool error = false;

	if (bits[0].rfind("challenge", 0) == 0)
	{
		try
		{
			next = doChallenge(lines, line);
		}
		catch (int e)
		{
			std::cout << "Unrecoverable error" << " ";
			std::cout << _filename << " : " << next << std::endl;
			exit(0);
		}
		catch (std::string err)
		{
			if (!_forgiving)
			{
				std::cout << "Unrecoverable error:" << std::endl;
				std::cout << _filename << " : " << next << std::endl;
				std::cout << err << std::endl;
				error = true;
				exit(0);
			}
		}
	}
	
	return next;
}

size_t Loader::doStrain(std::vector<std::string > &lines, size_t line)
{
	std::vector<std::string> bits = split(lines[line], ',');

	if (bits.size() <= 1)
	{
		std::cout << "Line " << line << ": STRAIN label found but "
		"no strain name given." << std::endl;
		throw -1;
	}
	
	std::string name = bits[1];
	
	bool check = false;
	if (_name2Strain.count(name) > 0)
	{
		std::cout << "WARNING: strain name " << name << " on line " << line <<
		" already exists. Will check that this is identical." << std::endl;
		check = true;
	}

	Strain *strain = new Strain(name);
	strain->setLoader(this);
	
	size_t next = line + 1;
	std::string mut_list;
	std::vector<std::string> sera_list;
	
	while (next < lines.size())
	{
		bits = split(lines[next], ',');
		std::string mName;
		std::string sName;

		if (bits.size() == 0 || (bits[0].length() > 0 && bits[0][0] == '#'))
		{
			next++;
			continue;
		}

		bits.resize(2);
		trim(bits[0]);
		trim(bits[1]);
		
		if (bits[0] == "strain" || bits[0] == "challenge")
		{
			break;
		}

		mName = bits[0];
		sName = bits[1];
		
		if (mName.length() > 0)
		{
			mut_list += mName + ",";
		}

		if (sName.length() > 0)
		{
			sera_list.push_back(sName);
		}
		
		next++;
	}

	strain->setList(mut_list, _muts);

	if (!check)
	{
		std::cout << strain->name() << " added" << std::endl;
		_name2Strain[strain->name()] = strain;
		_strains.push_back(strain);

		addNewSera(strain, sera_list);
	}
	else
	{
		Strain *existing = _name2Strain[name];
		if (!strain->sameAsStrain(existing))
		{
			std::cout << strain->name() << ": " << strain->summary();
			std::cout << std::endl;
			std::cout << existing->name() << ": " << existing->summary();
			std::cout << std::endl;
			std::cout << "Duplicate but non-matching definitions: " <<
			strain->name() << " and " << existing->name() << std::endl;
			throw -1;
		}
	}
	
	return next;
}

void Loader::addNewSera(Strain *str, std::vector<std::string> &list)
{
	for (size_t i = 0; i < list.size(); i++)
	{
		trim(list[i]);
		if (list[i].length() == 0)
		{
			continue;
		}

		for (size_t j = 0; j < _sera.size(); j++)
		{
			if (_sera[j]->name() == list[i])
			{
				std::cout << "Serum (" << list[i] << ") has already been "
				"associated with a strain (" << _sera[j]->strain()->name() 
				<< ")!" << std::endl;
				throw(-1);
			}
		}

		Serum *serum = new Serum(list[i], str);
		str->addChild(serum);
		_sera.push_back(serum);
		_name2Serum[serum->name()] = serum;
	}
}

void Loader::silenceMutations(std::string str)
{
	trim(str);
	std::vector<std::string> mutations = split(str, ',');
	
	for (size_t i = 0; i < mutations.size(); i++)
	{
		for (size_t j = 0; j < _muts.size(); j++)
		{
			if (_muts[j]->str() == mutations[i])
			{
				_muts[j]->silence(true);
			}
		}
	}
}

bool Loader::isAcceptable(int res)
{
	if (res == -1 || _acceptResidues.size() == 0)
	{
		return true;
	}
	
	for (size_t i = 0; i < _acceptResidues.size(); i++)
	{
		if (_acceptResidues[i] == res)
		{
			return true;
		}
	}
	
	return false;
}
