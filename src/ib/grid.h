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

#ifndef IB_GRID_H
#define IB_GRID_H

#include <cmath>
#include <type_traits>

#include "dt/point.h"
#include "dt/vector.h"

#define pi 3.141592653589793238462643383279502884L

namespace ib {

	template<class Geometry, unsigned int Size>
	class grid : public dt::vector<Size> {
		public:
			static constexpr unsigned int n = Geometry::n;
			static constexpr unsigned int nx = Geometry::nx;
			static constexpr unsigned int ny = Geometry::ny;
			static constexpr unsigned int nz = Geometry::nz;

			typedef Geometry geometry_type;
			typedef dt::point point_type;
			typedef typename point_type::value_type coordinate_type;
			using typename dt::vector<Size>::value_type;
		protected:
			virtual point_type coordinates(unsigned int) = 0;
			virtual std::vector<point_type> search_nearby(point_type) = 0;
			virtual double volume(unsigned int) const = 0;
			virtual double delta_function(coordinate_type r) {
				if (abs(n * r) > 2.)
					return 0;
				return 0.25 * n * (1 + cos(pi * r * n / 2));
			}
		public:
			template<unsigned int R>
				void interpolate(grid<Geometry, R>& source);
			grid(const grid& g) : dt::vector<Size>(g) {}
			grid(dt::vector<Size>& v) : dt::vector<Size>(v) {}
			grid(value_type v = 0.0) : dt::vector<Size>(v) {}
			grid(std::initializer_list<value_type> v) : dt::vector<Size>(v) {};
			using dt::vector<Size>::operator[];

		template<class, unsigned int> friend class grid;
	};

	template<class Geometry, unsigned int Size> template<unsigned int Size2>
	void
	grid<Geometry, Size>::interpolate(grid<Geometry, Size2>& source)
	{
		for (unsigned int i = 0; i < Size; ++i) {
			auto coords = coordinates(i);
			value_type value = value_type();
			for (auto& p: source.search_nearby(coords)) {
				auto index = p.get_index();
				auto weight = delta_function(p[0] - coords[0])
				            * delta_function(p[1] - coords[1])
				            * delta_function(p[2] - coords[2]);
				value += weight * source.volume(index) * source[index];
			}
			operator[](i) = value;
		}
	}
}

#endif
