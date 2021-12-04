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

#include "Tent.h"
#include "Plotter.h"
#include "Workbench.h"
#include "Scatter.h"
#include "Strain.h"
#include "Mutation.h"
#include <c4xsrc/GLAxis.h>
#include <QTreeWidget>
#include <h3dsrc/Icosahedron.h>

Plotter::Plotter(QWidget *parent) : SlipGL(parent)
{
	_tent = NULL;
	_r = 1; _g = 1; _b = 1;
	setFocus(Qt::MouseFocusReason);
	setFocusPolicy(Qt::StrongFocus);
	_controlPressed = false;
	_shiftPressed = false;
	
	if (false)
	{
		GLAxis *x = new GLAxis(make_vec3(1, 0, 0));
		x->addText("x");
		addObject(x);
	}

	if (false)
	{
		GLAxis *y = new GLAxis(make_vec3(0, 1, 0));
		y->addText("y");
		addObject(y);
	}
	
	Icosahedron *ico = new Icosahedron();
	ico->triangulate();
	ico->triangulate();
	ico->resize(0.02);
	addObject(ico);
	
	int dims = Workbench::dimensions();
	
	if (dims > 2)
	{
		GLAxis *z = new GLAxis(make_vec3(0, 0, 1));
		z->addText("z");
		addObject(z);
	}
	
	_scatter = new Scatter();
	_depth = false;
	zFar = 100;
	addObject(_scatter);
}

void Plotter::initializeGL()
{
	SlipGL::initializeGL();
	glDisable(GL_DEPTH_TEST);
}

void Plotter::setDepth(bool on)
{
	initializeGL();

	_depth = on;

	if (on)
	{
		std::cout << "Enabling" << std::endl;
		glEnable(GL_DEPTH_TEST);
		glDepthMask(GL_TRUE);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
	else
	{
		std::cout << "Disabling" << std::endl;
		glDisable(GL_DEPTH_TEST);
	}

}

void Plotter::setTree(QTreeWidget *w)
{
	_tree = w;
	connect(_tree, &QTreeWidget::itemSelectionChanged,
	        this, &Plotter::replotSelection);

	_tree->setSelectionMode(QAbstractItemView::ExtendedSelection);
}

void Plotter::replot(const std::vector<Mutation *> &muts)
{
	_scatter->populateFromMutations(muts);

	for (size_t i = 0; i < muts.size(); i++)
	{
		if (muts[i]->silenced())
		{
			continue;
		}

		_tree->addTopLevelItem(muts[i]);
	}
}

void Plotter::replot(const std::vector<Strain *> &strains)
{
	_scatter->populateFromStrains(strains);

	for (size_t i = 0; i < strains.size(); i++)
	{
		_tree->addTopLevelItem(strains[i]);
	}
}

void Plotter::replotSelection()
{
	std::vector<Mutation *> list;
	std::vector<Strain *> strains;
	QList<QTreeWidgetItem *> items = _tree->selectedItems();

	for (size_t i = 0; i < items.size(); i++)
	{
		Mutation *m = dynamic_cast<Mutation *>(items[i]);
		if (m != NULL)
		{
			list.push_back(m);
		}

		Strain *s = dynamic_cast<Strain *>(items[i]);
		if (s != NULL)
		{
			strains.push_back(s);
		}
	}
	
	if (items.size() == 0)
	{
		for (size_t j = 0; j < _tree->topLevelItemCount(); j++)
		{
			Mutation *m = dynamic_cast<Mutation *>(_tree->topLevelItem(j));
			Strain *s = dynamic_cast<Strain *>(_tree->topLevelItem(j));
			
			if (m != NULL)
			{
				list.push_back(m);
			}
			
			if (s != NULL)
			{
				strains.push_back(s);
			}
		}
	}

	_scatter->populateFromMutations(list);
	if (strains.size())
	{
		_scatter->populateFromStrains(strains);
	}
}

void Plotter::setVisibleText(bool vis)
{
	_scatter->setShowText(vis);
	
	for (size_t i = 0; i < 16; i++)
	{
		std::cout << _model.vals[i] << " ";
	}
	
	for (size_t i = 0; i < 3; i++)
	{
		std::cout << *(&_centre.x + i) << " ";
	}
	
	std::cout << std::endl;
}

void Plotter::setTent(Tent *t)
{
	if (_tent != NULL)
	{
		removeObject(_tent);
		delete _tent;
		_tent = NULL;
	}

	_tent = t;
	
	if (_tent)
	{
		addObject(_tent);
	}

	initializeGL();
	setDepth(_depth);
}

void Plotter::setShowsText(bool show)
{
	_scatter->setShowText(show);
}

void Plotter::setHeatMode(bool heat)
{
	Tent::setHeatMode(heat);
	
	if (_tent != NULL)
	{
		_tent->useMode();
	}
	
	if (!heat && _tent)
	{
		for (size_t i = 0; i < _objects.size(); i++)
		{
			_objects[i]->setDisabled(true);
		}
		
		_tent->setDisabled(false);
	}
	else
	{
		for (size_t i = 0; i < _objects.size(); i++)
		{
			_objects[i]->setDisabled(false);
		}
	}
	
	initializeGL();
	setDepth(true);
}
