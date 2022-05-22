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

#include "Display.h"
#include <QVBoxLayout>
#include <QMenu>
#include <QMenuBar>
#include <QCheckBox>
#include <QPushButton>
#include <QMouseEvent>
#include <QTabWidget>
#include <h3dsrc/Dialogue.h>
#include <hcsrc/FileReader.h>
#include <iostream>
#include "SerumView.h"
#include "Plotter.h"
#include <QTreeWidget>

#define deg2rad(a) (M_PI * (a) / 180)

Display::Display(int argc, char *argv[])
{
	QWidget *central = new QWidget(this);
	QVBoxLayout *vbox = new QVBoxLayout(NULL);
	central->setLayout(vbox);

	_tabs = new QTabWidget(central);
	vbox->addWidget(_tabs);
	
	SerumView *view = new SerumView(NULL);
	_view = view;
	_tabs->addTab(view, "Serology matrix");

	Plotter *mutPlot, *strainPlot;

	{
		QWidget *pp = new QWidget(NULL);
		QHBoxLayout *hbox = new QHBoxLayout(NULL);
		pp->setLayout(hbox);

		QTreeWidget *w = new QTreeWidget(NULL);
		hbox->addWidget(w);
		w->setMaximumSize(200, 10000);

		mutPlot = new Plotter(NULL);
		hbox->addWidget(mutPlot);
		mutPlot->setTree(w);
		mutPlot->setMinimumSize(1000, 900);

		_tabs->addTab(pp, "Mutation plot");
	}

	{
		QWidget *pp = new QWidget(NULL);
		QHBoxLayout *hbox = new QHBoxLayout(NULL);
		pp->setLayout(hbox);

		QVBoxLayout *vbox = new QVBoxLayout(NULL);
		QTreeWidget *w = new QTreeWidget(NULL);
		vbox->addWidget(w);
		w->setMaximumSize(200, 10000);
		
		QPushButton *p = new QPushButton("Reset view", NULL);
		connect(p, &QPushButton::clicked, _view, &SerumView::resetView);
		vbox->addWidget(p);
		
		QCheckBox *ch = new QCheckBox("Show text", NULL);
		ch->setChecked(true);
		ch->setObjectName("show_text");
		connect(ch, &QCheckBox::stateChanged, this, &Display::toggleText);
		vbox->addWidget(ch);

		hbox->addLayout(vbox);

		strainPlot = new Plotter(NULL);
		hbox->addWidget(strainPlot);
		strainPlot->setTree(w);
		strainPlot->setMinimumSize(1000, 900);
		w->setContextMenuPolicy(Qt::CustomContextMenu);
		connect(w, &QTreeWidget::customContextMenuRequested,
		        this, &Display::rightClickStrain);
		_tabs->addTab(pp, "Strain plot");
	}

	view->setPlotters(mutPlot, strainPlot);

	view->setCommandLineArgs(argc, argv);

	_tabs->show();
	
	setCentralWidget(central);
	makeMenu();
}

void Display::rightClickStrain()
{

}

void Display::makeMenu()
{
	QMenu *render = menuBar()->addMenu(tr("&Render"));
	QAction *act = render->addAction(tr("Rotating movie"));
	connect(act, &QAction::triggered, this, &Display::rotatingMovie);

	act = render->addAction(tr("Rotate 90 around Z"));
	connect(act, &QAction::triggered, this, &Display::rotateZ);
	act = render->addAction(tr("Rotate 90 around X"));
	connect(act, &QAction::triggered, this, &Display::rotateX);
}

void Display::toggleText()
{
	QCheckBox *cb = findChild<QCheckBox *>("show_text");
	bool checked = cb->checkState();

	_view->strainPlot()->setVisibleText(checked);
}

void Display::rotateX()
{
	_view->strainPlot()->rotate(deg2rad(-90), 0, 0);
}

void Display::rotateZ()
{
	_view->strainPlot()->rotate(0, deg2rad(-90), 0);
}

void Display::rotatingMovie()
{
	std::string folder = openDialogue(_view, "Folder for movie images", 
	                                  "", true, true);

	if (folder == "")
	{
		return;
	}

	std::string pattern = folder + "/fr*.png";
	std::vector<std::string> files = glob(pattern);

	for (size_t i = 0; i < files.size(); i++)
	{
		remove(files[i].c_str());
	}

	FileReader::setOutputDirectory(folder);

	for (size_t i = 0; i < 360; i++)
	{
		std::string number = i_to_str(i);
		std::string zeros = std::string(5 - number.length(), '0');
		std::string filename = "fr_" + zeros + number + ".png";
		std::string path = FileReader::addOutputDirectory(filename);
		std::cout << path << std::endl;

		_view->strainPlot()->rotate(0, deg2rad(-1), 0);
		_view->strainPlot()->saveImage(path);
	}
}
