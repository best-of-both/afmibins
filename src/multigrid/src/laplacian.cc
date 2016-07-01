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
#include "types/matrix.h"
#include "laplacian.h"

namespace multigrid {

	using types::size_type;
	using types::index_type;

	laplacian::laplacian(const types::geometry& geometry) :
		types::sparse_square_matrix(geometry.nx * geometry.ny, 5),
		m_geometry(geometry)
	{
		const size_type &nx = m_geometry.nx, &ny = m_geometry.ny,
		                &n = m_geometry.n;
		for (index_type row = 0u; row < nx * ny; ++row) {
			index_type x = row % nx, y = row / nx;
			const double offdiagonal = n * n;
			double diagonal = -4.0 * offdiagonal,
			       left = offdiagonal,
			       right = offdiagonal,
			       top = offdiagonal,
			       bottom = offdiagonal;

			if (nx == 1) {
				diagonal += left + right;
				left = right = 0;
			} else if (nx == 2) {
				if (x == 0) {
					right += left;
					left = 0;
				} else {
					left += right;
					right = 0;
				}
			}

			if (y == 0) {
				diagonal += bottom;
				bottom = 0;
			}
			if (y == ny - 1) {
				diagonal += top;
				top = 0;
			}

			if (y > 0)       push_value(row, row - nx, bottom);
			if (x == nx - 1) push_value(row, nx * y, right);
			if (x > 0)       push_value(row, row - 1, left);
			                 push_value(row, row, diagonal);
			if (x < nx - 1)  push_value(row, row + 1, right);
			if (x == 0)      push_value(row, nx - 1 + nx * y, left);
			if (y < ny - 1)  push_value(row, row + nx, top);
		}
	}

}
