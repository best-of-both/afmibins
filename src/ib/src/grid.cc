/*************************************************************************\
 *                                                                       *
 * Immersed Boundary Incompressible Navier-Stokes solver                 *
 *                                                                       *
 * Copyright (C) 2016  Andrew Kassen <atkassen@gmail.com>                *
 *                                                                       *
 * This program is free software: you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or     *
 * (at your option) any later version.                                   *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. *
 *                                                                       *
\*************************************************************************/

#include "types/typedefs.h"
#include "types/point.h"
#include "grid.h"

namespace ib {

	using types::index_type;

	void
	grid::interpolate(grid& source)
	{
		for (index_type i = 0; i < size(); ++i) {
			types::point coords = coordinates(i);
			double value = 0;
			for (types::point& p: source.search_nearby(coords)) {
				index_type index = p.index;
				types::point diff = p - coords;
				double weight = delta(diff.x) * delta(diff.y);
				value += weight * source.volume(index) * source[index];
			}
			types::vector::operator[](i) = value;
		}
	}

}
