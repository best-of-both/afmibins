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

#ifndef IB_MAC_GRID_H
#define IB_MAC_GRID_H

#include <vector>
#include "eulerian_grid.h"

namespace ib {

	namespace mac { enum component { x, y, z }; }

	template<class Geometry, mac::component Component>
	class mac_grid : public eulerian_grid<Geometry, Geometry::nx * Geometry::ny * Geometry::nz> {
		private:
			static constexpr unsigned int Size = Geometry::nx * Geometry::ny * Geometry::nz;
		public:
			using eulerian_grid<Geometry, Size>::n;
			using typename eulerian_grid<Geometry, Size>::point_type;
			using typename eulerian_grid<Geometry, Size>::value_type;
		protected:
			virtual std::vector<point_type> search_nearby(point_type);
		public:
			virtual point_type coordinates(unsigned int);
			mac_grid(value_type v = value_type()) : eulerian_grid<Geometry, Size>(v) {}
			mac_grid(std::initializer_list<value_type> v)
				: eulerian_grid<Geometry, Size>(v) {}
	};

	template<class Geometry, mac::component Component>
	auto
	mac_grid<Geometry, Component>::coordinates(unsigned int i) ->
		point_type
	{
		point_type point = eulerian_grid<Geometry, Size>::coordinates(i);
		point[(Component+1) % 3] += 0.5 / n;
		point[(Component+2) % 3] += 0.5 / n;
		return point;
	}

	template<class Geometry, mac::component Component>
	auto
	mac_grid<Geometry, Component>::search_nearby(point_type point) ->
		std::vector<point_type>
	{
		point[(Component+1) % 3] -= 0.5 / n;
		point[(Component+2) % 3] -= 0.5 / n;
		auto results = eulerian_grid<Geometry, Size>::search_nearby(point);
		for (auto& p: results) {
			p[(Component+1) % 3] += 0.5 / n;
			p[(Component+2) % 3] += 0.5 / n;
		}
		return results;
	}

	template<class Geometry>
	class mac_grid<Geometry, mac::y> : public eulerian_grid<Geometry, Geometry::nx * (Geometry::ny - 1) * Geometry::nz> {
		private:
			static constexpr unsigned int Size = Geometry::nx * (Geometry::ny - 1) * Geometry::nz;
		public:
			using eulerian_grid<Geometry, Size>::nx;
			using eulerian_grid<Geometry, Size>::ny;
			using eulerian_grid<Geometry, Size>::n;
			using eulerian_grid<Geometry, Size>::width;
			using eulerian_grid<Geometry, Size>::height;
			using eulerian_grid<Geometry, Size>::depth;
			using typename eulerian_grid<Geometry, Size>::point_type;
			using typename eulerian_grid<Geometry, Size>::value_type;
		protected:
			virtual std::vector<point_type> search_nearby(point_type);
		public:
			virtual point_type coordinates(unsigned int);
			mac_grid(value_type v = value_type()) : eulerian_grid<Geometry, Size>(v) {}
			mac_grid(std::initializer_list<value_type> v)
				: eulerian_grid<Geometry, Size>(v) {}
	};

	template<class Geometry>
	auto
	mac_grid<Geometry, mac::y>::coordinates(unsigned int i) ->
		point_type
	{
		return point_type(((double) (i % nx) + 0.5) / n,
		                  ((double) (i / nx %  (ny - 1)) + 1.0) / n,
		                  ((double) (i / (nx * (ny - 1))) + 0.5) / n, i);
	}

	template<class Geometry>
	auto
	mac_grid<Geometry, mac::y>::search_nearby(point_type coords) ->
		std::vector<point_type>
	{
		std::vector<point_type> results;

		coords[0] = std::fmod(std::fmod(coords[0], width) + width, width);
		coords[2] = std::fmod(std::fmod(coords[2], depth) + depth, depth);

		const auto fx = std::floor(coords[0] * n - 0.5),
				   fy = std::floor(coords[1] * n - 1.0),
				   fz = std::floor(coords[2] * n - 0.5);
		for (auto i = 0; i < 64; ++i) {
			const auto dix = i % 4,
			           diy = (i / 4) % 4,
			           diz = (i / 16) % 4;
			const auto ix = ((std::size_t) fx + dix + Geometry::nx - 1) % Geometry::nx,
			           iy = (std::size_t) fy + diy,
			           iz = ((std::size_t) fz + diz + Geometry::nz - 1) % Geometry::nz;
			const auto x = (fx + dix - 1.0) / n,
			           y = (fy + diy) / n,
			           z = (fz + diz - 1.0) / n;
			if (y < 0 || y > height)
				continue;
			results.push_back(point_type(x + 0.5 / n, y, z + 0.5 / n, ix + nx * (iy - 1 + (ny - 1) * iz)));
		}

		return results;
	}

}

#endif
