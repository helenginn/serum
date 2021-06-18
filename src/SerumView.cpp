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

#include "SerumView.h"
#include "MyDictator.h"
#include "Loader.h"
#include <iostream>
#include <QVBoxLayout>
#include <QMouseEvent>
#include <QLabel>
#include <c4xsrc/MatrixView.h>
//#include <QChart>
//#include <QChartView>

SerumView::SerumView(QWidget *p) : QMainWindow(p)
{
	_worker = NULL;
	_loader = new Loader();
	
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
	for (int i = 1; i < argc; i++)
	{
		std::string str = argv[i];
		_args.push_back(str);
	}

	_dictator = new MyDictator(this);
	_dictator->setArgs(_args);
	_dictator->run();
}

void SerumView::loadDefinitions(std::string filename)
{
	_loader->load(filename);
	updateView();
}

void SerumView::resultVectors(std::string filename)
{
	_loader->writeResultVectors(filename);
	
}

void SerumView::refine()
{
	bool ready = prepareWorkForObject(_loader);
	
	if (!ready)
	{
		return;
	}

	connect(_loader, SIGNAL(resultReady()), this, SLOT(handleResult()));

	_loader->refineLoop(25);
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
	_loader->populateRaw(&raw, 0);
	_loader->populateNames(&names);
	
	double w = _loader->serumCount();
	double h = _loader->strainCount();

	MatrixView *mv = new MatrixView(NULL, w * 5, h * 5);
	mv->populate(w, h, raw);
	_dataLabel->setPixmap(QPixmap::fromImage(*mv));

	_loader->populateRaw(&raw, 1);
	mv->populate(w, h, raw);
	_modelLabel->setPixmap(QPixmap::fromImage(*mv));

	_loader->populateRaw(&raw, -1);
	mv->populate(w, h, raw);
	_errorLabel->setPixmap(QPixmap::fromImage(*mv));
	_mv = mv;
	_mv->setNames(names);
	
	free(raw);
}

void SerumView::writeOut(std::string filename, int type)
{	
	_loader->writeOut(filename, type);
}

void SerumView::antigenicity(std::string filename)
{	
	_loader->antigenicity(filename);
}

void SerumView::setScale(double scale)
{
	_loader->setScale(scale);
	
}

void SerumView::handleResult()
{
	_worker->quit();
	updateView();
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
	bool ready = prepareWorkForObject(_loader);
	
	if (!ready)
	{
		return;
	}

	connect(_loader, SIGNAL(resultReady()), this, SLOT(handleResult()));

	for (size_t i = 0; i < 30; i++)
	{
		_loader->refineLoop(25);
		if (i % 3 == 0)
		{
			_loader->writeResultVectors("vectors.csv");
			_loader->antigenicity("antigen.csv");
		}
	}

	_loader->writeResultVectors("vectors.csv");
}
