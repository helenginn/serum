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

#include "MyDictator.h"
#include "SerumView.h"
#include "Loader.h"

MyDictator::MyDictator(SerumView *s) : Dictator()
{
	_serum = s;
}

bool MyDictator::processRequest(std::string first, std::string last)
{
	if (first == "quit")
	{
		exit(0);
	}
	else if (first == "load")
	{
		_serum->loadDefinitions(last);
	}
	else if (first == "remove-unused")
	{
		_serum->removeUnused();
	}
	else if (first == "write-errors")
	{
		_serum->writeOut(last, -1);
	}
	else if (first == "write-model")
	{
		_serum->writeOut(last, 1);
	}
	else if (first == "write-data")
	{
		_serum->writeOut(last, 0);
	}
	else if (first == "result-vectors")
	{
		_serum->resultVectors(last);
	}
	else if (first == "write-antigenicity")
	{
		_serum->antigenicity(last);
	}
	else if (first == "per-residue")
	{
		_serum->loader()->perResidueAntigenicity(last);
	}
	else if (first == "scale")
	{
		_serum->setScale(atof(last.c_str()));
	}
	else if (first == "refine")
	{
		_serum->refine();
	}
	else if (first == "average-model")
	{
		_serum->loader()->setModelToAverage();
	}
	else if (first == "refine-offsets")
	{
		_serum->loader()->setRefineOffset(true);
	}
	else if (first == "refine-strengths")
	{
		_serum->loader()->setRefineStrength(true);
	}
	else if (first == "dimension")
	{
		_serum->loader()->setDimension(atoi(last.c_str()));
	}
	else if (first == "run")
	{
		_serum->run();
	}

	return true;
}


