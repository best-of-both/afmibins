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

#include <cmath>

#include "point.h"
#include "geometry.h"

namespace types {

	point
	geometry::fit_to_box(point p) const
	{
		using std::fmod;
		p.x = fmod(fmod(p.x, width) + width, width);
		p.z = fmod(fmod(p.z, depth) + depth, depth);

		return p;
	}

	bool
	geometry::fits_to_box(point p) const
	{
		return p.y > 0 && p.y < height;
	}

}
