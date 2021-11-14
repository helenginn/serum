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
#include "Projection.h"
#include <h3dsrc/Icosahedron.h>
#include "Workbench.h"
#include "SerumView.h"
#include "Strain.h"
#include "MyDictator.h"
#include "DataPlot.h"
#include "Loader.h"
#include <iostream>
#include <QVBoxLayout>
#include <QMouseEvent>
#include <QLabel>
#include <c4xsrc/MatrixView.h>
#include <hcsrc/FileReader.h>
#include "Plotter.h"

SerumView::SerumView(QWidget *p) : QMainWindow(p)
{
	_dataPlot = NULL;
	_mutPlot = NULL;
	_strainPlot = NULL;
	_worker = NULL;
	_bench = new Workbench();
	_loader = new Loader(_bench);
	
	QWidget *central = new QWidget(this);
	QVBoxLayout *vbox = new QVBoxLayout(NULL);
	central->setLayout(vbox);

	{
		QLabel *l = new QLabel("Data: ", this);
		vbox->addWidget(l);
	}

	_dataPlot = new DataPlot(this);
	_dataPlot->setType(0);
	_dataPlot->setBench(_bench);
	vbox->addWidget(_dataPlot);

	{
		QLabel *l = new QLabel("Model: ", this);
		vbox->addWidget(l);
	}

	_modelPlot = new DataPlot(this);
	_modelPlot->setType(1);
	_modelPlot->setBench(_bench);
	vbox->addWidget(_modelPlot);

	{
		QLabel *l = new QLabel("Error: ", this);
		vbox->addWidget(l);
	}

	_errorPlot = new DataPlot(this);
	_errorPlot->setType(-1);
	_errorPlot->setBench(_bench);
	vbox->addWidget(_errorPlot);
	
	setCentralWidget(central);
}

void SerumView::setCommandLineArgs(int argc, char *argv[])
{
	show();

	for (int i = 1; i < argc; i++)
	{
		std::string str = argv[i];
		_args.push_back(str);
	}

	_dictator = new MyDictator(this);
	_dictator->setArgs(_args);
	_dictator->run();
}

void SerumView::loadSpreadsheet(std::string folder)
{
	_loader->parseSpreadsheet(folder);
	_bench->setup();
	updateView();
}

void SerumView::resultVectors(std::string filename)
{
	_bench->writeResultVectors(filename);
	
}

void SerumView::refine()
{
	bool ready = prepareWorkForObject(_bench);
	
	if (!ready)
	{
		return;
	}

	connect(_bench, SIGNAL(resultReady()), this, SLOT(handleResult()),
	        Qt::UniqueConnection);
	connect(this, SIGNAL(start()), _bench, SLOT(refineLoop()),
	        Qt::UniqueConnection);
	_worker->start();

	emit start();
//	_loader->refineLoop();
}

void SerumView::updateView()
{
	std::vector<Strain *> used;
	for (size_t i = 0; i < _bench->strainCount(); i++)
	{
		if (_bench->strain(i)->serumCount() > 0)
		{
			used.push_back(_bench->strain(i));
		}
	}

	_dataPlot->setStrains(used);
	_dataPlot->redraw();

	_modelPlot->setStrains(used);
	_modelPlot->redraw();

	_errorPlot->setStrains(used);
	_errorPlot->redraw();
	
	if (_mutPlot)
	{
		_mutPlot->replot(_bench->mutations());
	}
	
	if (_strainPlot)
	{
		_strainPlot->replot(_bench->strains());
	}
}

void SerumView::writeOut(std::string filename, int type)
{	
	_bench->writeOut(filename, type);
}

void SerumView::antigenicity(std::string filename)
{	
	_bench->antigenicity(filename);
}

void SerumView::setScale(double scale)
{
	_bench->setScale(scale);
	
}

void SerumView::handleResult()
{
	_worker->quit();
	_worker->wait();

	updateView();
	
	_dictator->incrementJob();
}

bool SerumView::prepareWorkForObject(QObject *object)
{
	if (object == NULL)
	{
		return false;
	}

	if (_worker && _worker->isRunning())
	{
		std::cout << "Waiting for worker to finish old job" << std::endl;
		_worker->wait();
	}
	
	if (!_worker)
	{
		_worker = new QThread();
	}

	object->moveToThread(_worker);

	return true;
}

void SerumView::run()
{
	bool ready = prepareWorkForObject(_bench);
	
	if (!ready)
	{
		return;
	}

	connect(_bench, SIGNAL(resultReady()), this, SLOT(handleResult()));

	for (size_t i = 0; i < 30; i++)
	{
		_bench->refineLoop();
		if (i % 3 == 0)
		{
			_bench->writeResultVectors("vectors.csv");
			_bench->antigenicity("antigen.csv");
		}
	}

	_bench->writeResultVectors("vectors.csv");
}

void SerumView::acceptMutations(std::string muts)
{
	std::vector<std::string> bits = split(muts, '-');
	int start = atoi(bits[0].c_str());

	if (bits.size() < 2)
	{
		_loader->addAcceptableResidue(start);
		return;
	}

	int end   = atoi(bits[1].c_str());
	
	for (int i = start; i < end; i++)
	{
		_loader->addAcceptableResidue(i);
	}
}

void SerumView::print(std::string prefix)
{
	_dataPlot->write(prefix);
	_modelPlot->write(prefix);
	_errorPlot->write(prefix);

}

void SerumView::tent(std::string strain)
{
	if (strain.length() == 0)
	{
		_strainPlot->setTent(NULL);
		return;
	}
	Tent *tent = new Tent(_bench);
	tent->makeBase(strain);
	std::cout << "Making tent" << std::endl;

	_strainPlot->setTent(tent);
}

void SerumView::resetView()
{
	changeView("");
}

void SerumView::changeView(std::string filename)
{
	if (filename.length() == 0)
	{
		filename = _view;
	}

	_view = filename;
	
	if (filename.length() == 0)
	{
		return;
	}

	Projection p(filename);
	p.applyToGL(_strainPlot);
	p.applyToGL(_mutPlot);
}

void SerumView::strainPhoto(std::string last)
{
	_strainPlot->saveImage(last);
}

void SerumView::mutationPhoto(std::string last)
{
	_mutPlot->saveImage(last);
}
