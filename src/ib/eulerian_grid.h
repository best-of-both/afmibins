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
#include "grid.h"

namespace ib {

	template<class Geometry, unsigned int Size = Geometry::nx * Geometry::ny * Geometry::nz>
	class eulerian_grid : public grid<Geometry, Size> {
		public:
			using grid<Geometry, Size>::n;
			using grid<Geometry, Size>::nx;
			using grid<Geometry, Size>::ny;
			using grid<Geometry, Size>::nz;
			using typename grid<Geometry, Size>::point_type;
			using typename grid<Geometry, Size>::value_type;
		protected:
			static constexpr unsigned int width = Geometry::width;
			static constexpr unsigned int height = Geometry::height;
			static constexpr unsigned int depth = Geometry::depth;

			virtual std::vector<point_type> search_nearby(point_type);
			virtual double volume(unsigned int) const { return 1. / (n * n * n); }
			virtual point_type coordinates(unsigned int);
		public:
			eulerian_grid(value_type v = value_type()) : grid<Geometry, Size>(v) {}
			eulerian_grid(grid<Geometry, Size>& v) : grid<Geometry, Size>(v) {}
			eulerian_grid(const eulerian_grid& v) : grid<Geometry, Size>(v) {}
			eulerian_grid(std::initializer_list<value_type> v) : grid<Geometry, Size>(v) {};
	};

	template<class Geometry, unsigned int Size>
	auto
	eulerian_grid<Geometry, Size>::coordinates(unsigned int index) ->
		point_type
	{
		return point_type(((double) (index % nx)) / n,
		                  ((double) (index % (nx * ny) / nx)) / n,
		                  ((double) (index / (nx * ny))) / n, index);
	}

	template<class Geometry, unsigned int Size>
	auto
	eulerian_grid<Geometry, Size>::search_nearby(point_type coords) ->
		std::vector<point_type>
	{
		std::vector<point_type> results;

		coords[0] = std::fmod(std::fmod(coords[0], width) + width, width);
		coords[2] = std::fmod(std::fmod(coords[2], depth) + depth, depth);

		for (auto i = 0; i < 64; ++i) {
			const auto dix = i % 4,
			           diy = (i / 4) % 4,
			           diz = (i / 16) % 4;
			const auto fx = std::floor(coords[0] * n),
			           fy = std::floor(coords[1] * n),
			           fz = std::floor(coords[2] * n);
			const auto ix = ((std::size_t) fx + dix + Geometry::nx - 1) % Geometry::nx,
			           iy = (std::size_t) fy + diy,
			           iz = ((std::size_t) fz + diz + Geometry::nz - 1) % Geometry::nz;
			const auto x = (fx + dix - 1.0) / n,
			           y = (fy + diy - 1.0) / n,
			           z = (fz + diz - 1.0) / n;
			if (y < 0 || y > height)
				continue;
			results.push_back(point_type(x, y, z, ix + nx * (iy - 1 + ny * iz)));
		}

		return results;
	}

}

#endif
