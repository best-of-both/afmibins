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

#include <vector>
#include <initializer_list>

#include "types/point.h"
#include "types/kdtree.h"
#include "grid.h"

namespace ib {

	class lagrangian_grid : public grid {
		protected:
			std::vector<types::point> m_positions;
			types::kdtree m_tree;
		public:
			virtual types::point coordinates(index_type index) const { return m_positions[index]; }
			virtual std::vector<types::point> search_nearby(types::point) const;

			lagrangian_grid(const types::geometry& geometry,
					std::vector<types::point> p, double v = 0.0) :
				grid(geometry, p.size(), v), m_positions(p), m_tree(p, geometry) {}
			lagrangian_grid(const lagrangian_grid& g) :
				grid(g), m_positions(g.m_positions), m_tree(g.m_positions, g.m_geometry) {}
			lagrangian_grid(const types::geometry& geometry,
					std::vector<types::point> p, std::initializer_list<double> v) :
				grid(geometry, v), m_positions(p), m_tree(p, geometry) {}
	};

}

#endif
