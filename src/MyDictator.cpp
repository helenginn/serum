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

#include <hcsrc/FileReader.h>
#include "MyDictator.h"
#include "SerumView.h"
#include "Settings.h"
#include "Workbench.h"
#include "Plotter.h"

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
	else if (first == "spreadsheet")
	{
		_serum->loadSpreadsheet(last);
	}
	else if (first == "free")
	{
		_serum->workbench()->makeStrainFree(last);
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
	else if (first == "cluster-mutations")
	{
		_serum->workbench()->clusterMutations(last);
	}
	else if (first == "write-antigenicity")
	{
		_serum->antigenicity(last);
	}
	else if (first == "per-residue")
	{
		_serum->workbench()->perResidueAntigenicity(last);
	}
	else if (first == "accept-residues")
	{
		_serum->acceptMutations(last);
	}
	else if (first == "scale")
	{
		_serum->setScale(atof(last.c_str()));
	}
	else if (first == "refine")
	{
		_serum->refine();
		return false;
	}
	else if (first == "refine-ease")
	{
		_serum->workbench()->setRefineEase(true);
	}
	else if (first == "refine-offsets")
	{
		_serum->workbench()->setRefineOffset(true);
	}
	else if (first == "refine-strengths")
	{
		_serum->workbench()->setRefineStrength(true);
	}
	else if (first == "refine-strain-strengths")
	{
		_serum->workbench()->setRefineStrainStrength(true);
	}
	else if (first == "refine-strain-offsets")
	{
		_serum->workbench()->setRefineStrainOffset(true);
	}
	else if (first == "refine-asymmetric")
	{
		_serum->workbench()->setRefineImportance(true);
	}
	else if (first == "dimension")
	{
		_serum->workbench()->setDimension(atoi(last.c_str()));
	}
	else if (first == "assign-to-pdb")
	{
		_serum->workbench()->assignToPDB(last);
	}
	else if (first == "tent-mode")
	{
		bool heat = false;
		if (last == "heat")
		{
			heat = true;
		}

		_serum->strainPlot()->setHeatMode(heat);
	}
	else if (first == "depth")
	{
		bool depth = true;
		if (last == "off")
		{
			depth = false;
		}

		_serum->strainPlot()->setDepth(depth);
		_serum->mutPlot()->setDepth(depth);
	}
	else if (first == "tent")
	{
		_serum->tent(last);
	}
	else if (first == "display-settings")
	{
		Settings s = _serum->workbench()->displaySettings(last.c_str());
		_serum->applySettings(s);
	}
	else if (first == "show-text")
	{
		_serum->strainPlot()->setShowsText(true);
	}
	else if (first == "hide-text")
	{
		_serum->strainPlot()->setShowsText(false);
	}
	else if (first == "strain-photo")
	{
		_serum->strainPhoto(last);
	}
	else if (first == "mutation-photo")
	{
		_serum->mutationPhoto(last);
	}
	else if (first == "select")
	{
		std::vector<Strain *> strains = _serum->workbench()->getStrains(last);
		_serum->strainPlot()->replot(strains);
	}
	else if (first == "view")
	{
		_serum->changeView(last);
	}
	else if (first == "run")
	{
		_serum->run();
	}
	else if (first == "random-seed")
	{
		srand(atoi(last.c_str()));
	}
	else if (first == "print")
	{
		_serum->print(last);
	}
	else if (first == "model-type")
	{
		to_lower(last);
		_properties[first] = last;
	}

	return true;
}


