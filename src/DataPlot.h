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

#ifndef __serum__DataPlot__
#define __serum__DataPlot__

#include <QGraphicsScene>
#include <QGraphicsView>
#include <vector>

class Strain;
class Workbench;
class Serum;
class MatrixView;

class DataPlot : public QGraphicsView
{
public:
	DataPlot(QWidget *parent);

	void setStrains(std::vector<Strain *> &strains);
	
	void setType(int type)
	{
		_type = type;
	}
	
	void setBench(Workbench *bench)
	{
		_bench = bench;
	}

	void redraw();
	void write(std::string prefix);
protected:
	virtual void resizeEvent(QResizeEvent *);
private:
	QGraphicsScene *_scene;
	
	std::vector<Strain *> _strains;
	std::vector<Serum *> _xsera, _ysera;
	int _type;

	Workbench *_bench;
};

#endif
