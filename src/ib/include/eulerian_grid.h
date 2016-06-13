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

#ifndef IB_EULERIAN_GRID_H
#define IB_EULERIAN_GRID_H

#include <vector>
#include <initializer_list>

#include "types/typedefs.h"
#include "types/point.h"
#include "grid.h"

namespace ib {

	using types::size_type;
	using types::index_type;

	class eulerian_grid : public grid {
		protected:
			using grid::m_geometry;

			virtual std::vector<types::point> search_nearby(types::point) const;
			virtual double volume(index_type) const;
			virtual types::point coordinates(index_type) const;
		public:
			eulerian_grid(const types::geometry& geometry, double v = 0.0) :
				grid(geometry, geometry.nx * geometry.ny * geometry.nz, v) {}
			eulerian_grid(const types::geometry& geometry, size_type size, double v = 0.0) :
				grid(geometry, size, v) {}
			eulerian_grid(const types::geometry& geometry, std::initializer_list<double> v) :
				grid(geometry, v) {}
	};

}

#endif
