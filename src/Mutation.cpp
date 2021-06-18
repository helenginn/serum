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

#include "Mutation.h"
#include "Strain.h"

Mutation::Mutation(std::string mut, int dims)
{
	_mut = mut;
	_dim = dims;
	_vec = new double[dims];
	for (size_t i = 0; i < dims; i++)
	{
		_vec[i] = 0;
	}
}

void Mutation::refresh(void *object)
{
	Mutation *mut = static_cast<Mutation *>(object);
	
	for (size_t i = 0; i < mut->_strains.size(); i++)
	{
		Strain *str = mut->_strains[i];
		
		str->needsRefresh();
	}
}

void Mutation::randomiseVector()
{
	for (size_t j = 0; j < _dim; j++)
	{
		_vec[j] = ((rand() / (double)RAND_MAX) - 0.5) / 2;
	}
}

std::vector<double> &Mutation::asVector()
{
	_vector.resize(_dim);
	for (size_t i = 0; i < _dim; i++)
	{
		_vector[i] = _vec[i];
	}

	return _vector;
}
