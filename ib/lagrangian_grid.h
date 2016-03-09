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

#ifndef IB_LAGRANGIAN_GRID_H
#define IB_LAGRANGIAN_GRID_H

#include <ostream>

#include "grid.h"
#include "dt/vector.h"
#include "dt/kdtree.h"

namespace ib {

	template<class Geometry, unsigned int Size>
	class lagrangian_grid : public grid<Geometry, Size> {
		public:
			using grid<Geometry, Size>::n;
			static constexpr unsigned int width = Geometry::width;
			static constexpr unsigned int height = Geometry::height;
			static constexpr unsigned int depth = Geometry::depth;

			using typename grid<Geometry, Size>::point_type;
			using typename grid<Geometry, Size>::value_type;
		private:
			std::vector<point_type> positions;
			dt::kdtree<Geometry> tree;
		public:
			virtual point_type coordinates(unsigned int index) { return positions[index]; }
			virtual std::vector<point_type> search_nearby(point_type);

			lagrangian_grid(std::vector<point_type> p, value_type v = value_type())
				: grid<Geometry, Size>(v), positions(p), tree(p) {}
			lagrangian_grid(const lagrangian_grid& g)
				: grid<Geometry, Size>(g), positions(g.positions), tree(g.positions) {}
			lagrangian_grid(std::vector<point_type> p, std::initializer_list<value_type> v)
				: grid<Geometry, Size>(v), positions(p), tree(p) {}
	};

	template<class Geometry, unsigned int np>
	auto
	lagrangian_grid<Geometry, np>::search_nearby(point_type coords) ->
		std::vector<point_type>
	{
		coords[0] = std::fmod(std::fmod(coords[0], width) + width, width);
		coords[2] = std::fmod(std::fmod(coords[2], depth) + depth, depth);

		return tree.search_nearby(coords);
	}

}

#endif
