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
#include "SerumView.h"
#include "MyDictator.h"
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
	_mutPlot = NULL;
	_strainPlot = NULL;
	_worker = NULL;
	_bench = new Workbench();
	_loader = new Loader(_bench);
	
	QWidget *central = new QWidget(this);
	QVBoxLayout *vbox = new QVBoxLayout(NULL);
	central->setLayout(vbox);

	{
		QHBoxLayout *hbox = new QHBoxLayout(NULL);
		QLabel *l = new QLabel("Data: ", this);
		_lWhat = new QLabel("", this);
		hbox->addWidget(l);
		hbox->addWidget(_lWhat);
		vbox->addLayout(hbox);
	}

	_dataLabel = new QLabel("", this);
	vbox->addWidget(_dataLabel);

	{
		QLabel *l = new QLabel("Model: ", this);
		vbox->addWidget(l);
	}

	_modelLabel = new QLabel("", this);
	vbox->addWidget(_modelLabel);

	{
		QLabel *l = new QLabel("Error: ", this);
		vbox->addWidget(l);
	}

	_errorLabel = new QLabel("", this);
	vbox->addWidget(_errorLabel);
	
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

void SerumView::mousePressEvent(QMouseEvent *e)
{
	std::string str;
	{
		QPoint p1 = _dataLabel->mapFromGlobal(e->globalPos());
		str = _mv->getName(p1.x(), p1.y());
	}

	if (str.length() == 0)
	{
		QPoint p1 = _modelLabel->mapFromGlobal(e->globalPos());
		str = _mv->getName(p1.x(), p1.y());
	}

	if (str.length() == 0)
	{
		QPoint p1 = _errorLabel->mapFromGlobal(e->globalPos());
		str = _mv->getName(p1.x(), p1.y());
	}
	
	_lWhat->setText(QString::fromStdString(str));
}

void SerumView::updateView()
{
	double **raw = NULL;
	char ***names = NULL;
	_bench->populateRaw(&raw, 0);
	_bench->populateNames(&names);
	
	double w = _bench->serumCount();
	double h = _bench->strainCount();

	MatrixView *mv = new MatrixView(NULL, w * 5, h * 5);
	mv->populate(w, h, raw);
	_dataLabel->setPixmap(QPixmap::fromImage(*mv));

	/*
	_bench->populateRaw(&raw, 1);
	mv->populate(w, h, raw);
	_modelLabel->setPixmap(QPixmap::fromImage(*mv));
	*/

	_bench->populateRaw(&raw, -1);
	mv->populate(w, h, raw);
	_errorLabel->setPixmap(QPixmap::fromImage(*mv));
	_mv = mv;
	_mv->setNames(names);
	
	free(raw);
	
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
