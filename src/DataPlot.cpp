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

#include "DataPlot.h"
#include "Strain.h"
#include "Serum.h"
#include "Workbench.h"
#include <iostream>
#include <QGraphicsItem>
#include <QPrinter>
#include <c4xsrc/MatrixView.h>

#define BOX_DIM 12
#define MARGIN (BOX_DIM * 2)

DataPlot::DataPlot(QWidget *parent) : QGraphicsView(parent)
{
	_scene = new QGraphicsScene(this);
	setScene(_scene);
	_type = 0;

	show();
}

void DataPlot::resizeEvent(QResizeEvent *)
{
	_scene->setSceneRect(0, 0, width(), height());
}

void DataPlot::setStrains(std::vector<Strain *> &strains)
{
	_strains = strains;
	_xsera.clear();
	_ysera.clear();
	
	for (size_t i = 0; i < _strains.size(); i++)
	{
		for (size_t j = 0; j < _strains[i]->childrenCount(); j++)
		{
			Serum *child = _strains[i]->child(j);
			if (std::find(_xsera.begin(), _xsera.end(), child) == _xsera.end())
			{
				_xsera.push_back(child);
			}
		}

		for (size_t j = 0; j < _strains[i]->serumCount(); j++)
		{
			Serum *child = _strains[i]->serum(j);
			if (std::find(_ysera.begin(), _ysera.end(), child) == _ysera.end())
			{
				_ysera.push_back(child);
			}
		}
	}
}

void DataPlot::redraw()
{
	qDeleteAll(_scene->items());
	_scene->clear();
	
	int left = 0;
	int bottom = 0;
	int font_height = 0;
	
	std::vector<QGraphicsTextItem *> labels;
	for (size_t n = 0; n < _strains.size(); n++)
	{
		QGraphicsTextItem *item;
		item = _scene->addText(QString::fromStdString(_strains[n]->name()));
		item->setPos(0, BOX_DIM * n);
		labels.push_back(item);
		
		int text_width = item->boundingRect().width();
		QFont f = item->font();
		QFontMetrics fm(f);
		font_height = fm.height();
		bottom = item->pos().y() + font_height;

		if (left < text_width)
		{
			left = text_width;
		}
	}

	for (size_t n = 0; n < _bench->strainCount(); n++)
	{
		if (_bench->strain(n)->childrenCount() == 0)
		{
			continue;
		}

		int h = _strains.size();
		int w = _bench->strain(n)->childrenCount();

		double **ptr = (double **)calloc(w, sizeof(double *));
		double *raw = (double *)calloc(w * h, sizeof(double));

		for (size_t i = 0; i < w; i++)
		{
			ptr[i] = &raw[i * h];
		}

		for (size_t i = 0; i < w; i++)
		{
			Serum *serum = _bench->strain(n)->child(i);

			for (size_t j = 0; j < h; j++)
			{
				Strain *str = _strains[j];
				double remove = serum->strain()->defaultOffset();

				if (!str->hasSerum(serum))
				{
					ptr[i][j] = NAN;
					continue;
				}

				double val = str->serumValue(serum) - remove;

				if (_type == 0)
				{
					ptr[i][j] = val;
				}

				double model = _bench->modelForPair(str, serum);
				model -= remove;

				if (_type == 1)
				{
					ptr[i][j] = model;
				}

				if (_type == -1)
				{
					ptr[i][j] = val - model;
				}

				ptr[i][j] *= 0.3;
				ptr[i][j] += 0.5;
			}
		}

		MatrixView *matrix = new MatrixView(NULL, 
		                                    w * BOX_DIM, h * BOX_DIM);
		matrix->populate(w, h, ptr);
		QPixmap pixmap = QPixmap::fromImage(*matrix);

		QGraphicsTextItem *item;
		item = _scene->addText(QString::fromStdString(_bench->strain(n)->name()));
		item->setPos(left, 0);
		int text_height = item->boundingRect().height();
		int text_width = item->boundingRect().width();

		QGraphicsPixmapItem *pix = _scene->addPixmap(pixmap);
		pix->setOffset(left, text_height);
		int shift = std::max(pixmap.width(), text_width);
		left += shift + MARGIN;
		
		int pixbot = pixmap.height() + text_height;
		int label_height_shift = pixbot - bottom - 2;
		
		for (size_t i = 0; i < labels.size(); i++)
		{
			labels[i]->setPos(labels[i]->pos().x(), labels[i]->pos().y()
			                  + label_height_shift);
		}
		labels.clear();

		free(ptr);
		free(raw);
		delete matrix;
	}
	
	_scene->setSceneRect(_scene->itemsBoundingRect());
	fitInView(_scene->sceneRect());
}

void DataPlot::write(std::string prefix)
{
	QString name = QString::fromStdString(prefix) + "_";
	if (_type == 0)
	{
		name += "data.pdf";
	}
	else if (_type == 1)
	{
		name += "model.pdf";
	}
	else if (_type == -1)
	{
		name += "error.pdf";
	}

	QPrinter printer(QPrinter::HighResolution);
	printer.setOutputFormat(QPrinter::PdfFormat);
	printer.setPageSize(QPageSize(_scene->sceneRect().size().toSize()));
	printer.setOutputFileName(name);

	QPainter painter(&printer);
	_scene->render(&painter);
}
