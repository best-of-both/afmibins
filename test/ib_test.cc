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

#include <vector>
#include <iostream>
#include <cmath>

#include "ib/lagrangian_grid.h"
#include "ib/eulerian_grid.h"
#include "ib/mac_grid.h"
#include "dt/geometry.h"
#include "dt/point.h"


template<class T>
class point : public ib::lagrangian_grid<T, 1u>
{
	private:
		static constexpr unsigned int n = T::n;
	public:
		using typename ib::lagrangian_grid<T, 1u>::point_type;
		double volume(unsigned int) const { return 1. / (n * n * n); }
		point(double x, double y, double z)
			: ib::lagrangian_grid<T, 1u>({dt::point(x, y, z, 0u)}) {}
		point& operator=(double v) { for (auto& e: *this) e = v; return *this; }
};

template<class T, unsigned int np>
class circle : public ib::lagrangian_grid<T, np>
{
	public:
		using typename ib::lagrangian_grid<T, np>::point_type;
	private:
		double radius;
		static std::vector<point_type> compute_points(double, double, double, double);
	public:
		double volume(unsigned int) const { return radius * 2 * pi / np; }
		circle(double x, double y, double z, double radius)
			: ib::lagrangian_grid<T, np>(compute_points(x, y, z, radius)), radius(radius) {}
};

template<class T, unsigned int np>
auto
circle<T, np>::compute_points(double x, double y, double z, double radius) ->
std::vector<point_type>
{
	std::vector<point_type> points;
	points.reserve(np);
	for (auto i = 0u; i < np; ++i) {
		auto angle = (2 * pi * i) / np;
		points.push_back(point_type(x + radius * cos(angle),
		                            y + radius * sin(angle), z, i));
	}
	return points;
}

int
main(void)
{
	typedef dt::geometry<5u, 1u, 1u, 1u> geometry;
	ib::mac_grid<geometry, ib::mac::y> egrid;
	point<geometry> lgrid(0.5, 0.5, 0.5);

	lgrid = 1.0;
	egrid.interpolate(lgrid);

	for (auto& e: egrid) {
		if (e != 0.0)
			e = 1;
	}

	lgrid = 0.0;
	lgrid.interpolate(egrid);
	std::cout << lgrid << std::endl;
	return 0;
}
