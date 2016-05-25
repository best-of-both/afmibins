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
#include <cmath>
#include <cassert>

#include "types/typedefs.h"
#include "types/geometry.h"
#include "types/point.h"
#include "eulerian_grid.h"
#include "lagrangian_grid.h"

using namespace ib;

class circle : public lagrangian_grid {
	private:
		const double m_radius;

		static std::vector<types::point> construct_circle(const types::point, const double);
	protected:
		using lagrangian_grid::m_positions;

		double volume(types::index_type) const { return 2 * M_PI * m_radius / m_positions.size(); }
	public:
		circle(const types::geometry& geometry, const types::point c, const double radius) :
			lagrangian_grid(geometry, construct_circle(c, radius), 1.0), m_radius(radius) {}
};

std::vector<types::point>
circle::construct_circle(const types::point center, const double radius)
{
	std::vector<types::point> results;
	results.reserve(360);
	for (unsigned int i = 0; i < 360; ++i) {
		const double alpha = 2 * M_PI / 360 * i;
		results.push_back({center.x + radius * cos(alpha),
				center.y + radius * sin(alpha), center.z, i});
	}
	return results;
}

int
main(void)
{
	const types::point center(0.5, 0.5, 0.5);
	const types::geometry geometry(1, 1, 1, 20);
	eulerian_grid egrid(geometry);
	circle lgrid(geometry, center, 0.25);
	
	egrid.interpolate(lgrid);

	double total = 0.0;
	for (auto v: egrid)
		total += v;

	assert(abs(total / (4000 * M_PI) - 1) < 1e-3);
}
