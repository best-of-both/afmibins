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

#ifndef IMMERSED_BOUNDARY_GRID_H
#define IMMERSED_BOUNDARY_GRID_H

#include <cmath>
#include <vector>

#include "types/typedefs.h"
#include "types/geometry.h"
#include "types/vector.h"
#include "types/point.h"

namespace ib {

	using types::size_type;
	using types::index_type;

	class grid : public types::vector {
		protected:
			const types::geometry& m_geometry;
			virtual types::point coordinates(index_type) const = 0;
			virtual std::vector<types::point> search_nearby(types::point) const = 0;
			virtual double volume(index_type) const = 0;
			virtual double delta(double r) const {
				const size_type &n = m_geometry.n;
				if (abs(n * r) > 2.)
					return 0;
				return 0.25 * n * (1 + cos(M_PI * r * n / 2));
			}
		public:
			void interpolate(grid& source);

			grid(const types::geometry& geometry, std::initializer_list<value_type> l) :
				types::vector(l), m_geometry(geometry) {}
			grid(const types::geometry& geometry, size_type size, double v = 0.0) :
				types::vector(size, v), m_geometry(geometry) {};
	};

}

#endif
