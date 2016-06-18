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
#include <cassert>
#include <cmath>

#include "types/geometry.h"
#include "point.h"
#include "kdtree.h"

using namespace types;

int
main(void)
{
	const types::geometry g(2, 2, 2, 2);
	std::vector<point> points = {
		{0, 0, 0, 0},
		{0, 0, 1, 1},
		{0, 1, 0, 2},
		{0, 1, 1, 3},
		{1, 0, 0, 4},
		{1, 0, 1, 5},
		{1, 1, 0, 6},
		{1, 1, 1, 7}
	};

	kdtree tree(points, g);
	point p(0.5, 0, 0, 0);

	for (auto& q: tree.search_nearby(p)) {
		point r = q - p;
		assert(abs(r.x) < 1);
	}

	return 0;
}
