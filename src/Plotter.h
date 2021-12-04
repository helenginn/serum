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
#ifndef __serum__Plotter__
#define __serum__Plotter__

#include <h3dsrc/SlipGL.h>

class QTreeWidget;
class Scatter;
class Mutation;
class Strain;
class Tent;

class Plotter : public SlipGL
{
Q_OBJECT
public:
	Plotter(QWidget *parent);
	
	void setTree(QTreeWidget *w);
	void setTent(Tent *t);
	
	void setShowsText(bool show);

	void replot(const std::vector<Mutation *> &muts);
	void replot(const std::vector<Strain *> &strains);
	
	void setVisibleText(bool vis);
	void setDepth(bool on);
	
	void setHeatMode(bool heat);
protected:
	virtual void initializeGL();
public slots:
	void replotSelection();
private:
	bool _depth;
	Scatter *_scatter;
	QTreeWidget *_tree;
	Tent *_tent;

};

#endif

