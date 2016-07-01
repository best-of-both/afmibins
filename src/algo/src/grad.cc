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
#include "types/vector.h"
#include "grad.h"

namespace algo {

	using types::size_type;
	using types::index_type;
	using types::vector;

	void
	grad(vector& px, vector& py, const types::geometry& geometry,
		const vector& p)
	{
		const size_type &nx = geometry.nx, &ny = geometry.ny,
		                &n = geometry.n;
		for (index_type i = 0u; i < nx * ny; ++i) {
			index_type x = i % nx, y = i / nx;

			px[i] = n * (p[i] - p[i + ((x + nx - 1) % nx - x)]);
			if (y > 0u)
				py[i - nx] = n * (p[i] - p[i - nx]);
		}
	}

}
