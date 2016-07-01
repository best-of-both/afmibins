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

#include "types/typedefs.h"
#include "types/geometry.h"
#include "types/vector.h"
#include "div.h"

namespace algo {

	using types::size_type;
	using types::index_type;
	using types::vector;

	void
	div(vector& d, const types::geometry& geometry,
			const vector& u, const vector& v,
			const vector& top, const vector& bottom)
	{
		const size_type &nx = geometry.nx, &ny = geometry.ny,
		                &n = geometry.n;
		for (index_type i = 0u; i < nx * ny; ++i) {
			index_type x = i % nx,
			           y = i / nx,
			           ixz = x;
			double dv;

			if (ny > 1) {
				if (y == 0u)
					dv = v[i] - bottom[ixz];
				else if (y == ny - 1)
					dv = top[ixz] - v[i - nx];
				else
					dv = v[i] - v[i - nx];
			}
			else
				dv = top[ixz] - bottom[ixz];
			d[i] = n * (
				(u[i + ((x + 1) % nx - x)] - u[i]) + dv
			);
		}
	}

}
