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

#ifndef __serum__SerumView__
#define __serum__SerumView__

#include <QMainWindow>
#include <QThread>

class MyDictator;
class Plotter;
class MatrixView;
class Loader;
class QThread;
class QLabel;
class QMouseEvent;
class Workbench;

class SerumView : public QMainWindow
{
Q_OBJECT
public:
	SerumView(QWidget *p);

	void setCommandLineArgs(int argc, char *argv[]);
	void loadSpreadsheet(std::string folder);
	void refine();
	void resultVectors(std::string filename);
	void writeOut(std::string filename, int type);
	void antigenicity(std::string filename);
	void setScale(double scale);
	void run();
	
	void acceptMutations(std::string muts);
	
	void setPlotters(Plotter *mut, Plotter *strain)
	{
		_mutPlot = mut;
		_strainPlot = strain;
	}
	
	Workbench *workbench()
	{
		return _bench;
	}
signals:
	void start();
protected:
	virtual void mousePressEvent(QMouseEvent *e);
private slots:
	void handleResult();
	void updateView();
private:
	bool prepareWorkForObject(QObject *object);
	MyDictator *_dictator;
	Loader *_loader;

	std::vector<std::string> _args;
	Workbench *_bench;
	MatrixView *_mv;
	QThread *_worker;
	QLabel *_dataLabel;
	QLabel *_modelLabel;
	QLabel *_errorLabel;
	QLabel *_lWhat;
	Plotter *_mutPlot;
	Plotter *_strainPlot;
};

#endif
