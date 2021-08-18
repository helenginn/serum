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

class Workbench;
class Mutation;
class Table;
class Strain;
class Serum;
class Any;

typedef std::complex<double> Complex;

class Loader : public QObject
{
Q_OBJECT
public:
	Loader(Workbench *b);

	void load(std::string filename);
	void parseSpreadsheet(std::string folder);
	
	void addAcceptableResidue(int i)
	{
		_acceptResidues.push_back(i);
	}
	
	bool isAcceptable(int i);
signals:
	void update();
private:
	void addNewSera(Strain *str, std::vector<std::string> &list);
	void loadCSV(std::string &filename);
	size_t doLine(std::vector<std::string > &lines, size_t line);
	size_t doStrain(std::vector<std::string > &lines, size_t line);
	size_t doChallenge(std::vector<std::string > &lines, size_t line);

	void silenceMutations(std::string str);

	Workbench *_bench;

	std::vector<Strain *> _strains;
	std::vector<Serum *> _sera;
	std::vector<Mutation *> _muts;
	std::vector<Table *> _tables;

	std::map<std::string, Strain *> _name2Strain;
	std::map<std::string, Serum *> _name2Serum;
	bool _forgiving;
	std::string _filename;
	std::string _folder;
	Table *_currentTable;
	std::vector<int> _acceptResidues;
};

#endif
