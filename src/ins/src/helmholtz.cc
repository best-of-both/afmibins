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
#include "helmholtz.h"

namespace ins {

	using types::size_type;
	using types::index_type;

	helmholtz_periodic::helmholtz_periodic(const types::geometry& geometry, double scale) :
		types::sparse_square_matrix(geometry.nx * geometry.ny * geometry.nz, 7),
		m_geometry(geometry), m_scale(scale)
	{
		const size_type &nx = geometry.nx, &ny = geometry.ny,
		                &nz = geometry.nz, &n = geometry.n;

		for (index_type row = 0; row < nx * ny * nz; ++row) {
			index_type x = row % nx, y = (row % (nx * ny)) / nx, z = row / (nx * ny);
			const double offdiagonal = scale * n * n;
			double diagonal = 1 - 6.0 * offdiagonal,
			       left = offdiagonal,
			       right = offdiagonal,
			       top = offdiagonal,
			       bottom = offdiagonal,
			       front = offdiagonal,
			       back = offdiagonal;

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

			if (ny == 1) {
				diagonal -= 3 * (top + bottom);
				top = bottom = 0;
			} else {
				if (y == 0) {
					diagonal += bottom;
					bottom = 0;
				}
				if (y == ny - 1) {
					diagonal += top;
					top = 0;
				}
			}

			if (nz == 1) {
				diagonal += back + front;
				back = front = 0;
			} else if (nz == 2) {
				if (z == 0) {
					front += back;
					back = 0;
				} else {
					back += front;
					front = 0;
				}
			}

			if (z == nz - 1) push_value(row, x + nx * y, front);
			if (z > 0)       push_value(row, row - nx * ny, back);
			if (y > 0)       push_value(row, row - nx, bottom);
			if (x == nx - 1) push_value(row, nx * (y + ny * z), right);
			if (x > 0)       push_value(row, row - 1, left);
			                 push_value(row, row, diagonal);
			if (x < nx - 1)  push_value(row, row + 1, right);
			if (x == 0)      push_value(row, nx - 1 + nx * (y + ny * z), left);
			if (y < ny - 1)  push_value(row, row + nx, top);
		 	if (z < nz - 1)  push_value(row, row + ny * ny, front);
			if (z == 0)      push_value(row, x + ny * (y + ny * (nz - 1)), back);
		}
	}

	types::vector
	helmholtz_periodic::boundary(types::vector& t, types::vector& b) const
	{
		const types::size_type &nx = m_geometry.nx, &ny = m_geometry.ny,
		                       &nz = m_geometry.nz, &n = m_geometry.n;
		types::vector r(nx * ny * nz);
		for (index_type x = 0; x < nx; ++x) {
			for (index_type z = 0; z < nz; ++z) {
				const index_type i = x + nx * z;
				if (ny > 1) {
					r[x + nx * ny * z] = 8. / 3. * m_scale * n * n * b[i];
					r[x + nx * (ny - 1 + ny * z)] = 8. / 3. * m_scale * n * n * t[i];
				}
				else {
					r[x + nx * ny * z] = 3 * m_scale * n * n * b[i];
					r[x + nx * (ny - 1 + ny * z)] = 3 * m_scale * n * n * t[i];
				}
			}
		}
		return r;
	}

	helmholtz_dirichlet::helmholtz_dirichlet(const types::geometry& geometry, double scale) :
		types::sparse_square_matrix(geometry.nx * (geometry.ny - 1) * geometry.nz, 7),
		m_geometry(geometry), m_scale(scale)
	{
		const size_type &nx = geometry.nx, &ny = geometry.ny,
		                &nz = geometry.nz, &n = geometry.n;

		for (index_type row = 0; row < nx * (ny - 1) * nz; ++row) {
			index_type x = row % nx, y = (row % (nx * (ny - 1))) / nx, z = row / (nx * (ny - 1));
			const double offdiagonal = m_scale * n * n;
			double diagonal = 1. - 6 * offdiagonal,
			       left = offdiagonal,
			       right = offdiagonal,
			       top = offdiagonal,
			       bottom = offdiagonal,
			       front = offdiagonal,
			       back = offdiagonal;

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

			if (y == 0) bottom = 0;
			if (y == ny - 2) top = 0;

			if (nz == 1) {
				diagonal += front + back;
				front = back = 0;
			} else if (nz == 2) {
				if (z == 0) {
					front += back;
					back = 0;
				} else {
					back += front;
					front = 0;
				}
			}

			if (z == nz - 1) push_value(row, x + nx * y, front);
			if (z > 0)       push_value(row, row - nx * ny, back);
			if (y > 0)       push_value(row, row - nx, bottom);
			if (x == nx - 1) push_value(row, nx * (y + ny * z), right);
			if (x > 0)       push_value(row, row - 1, left);
			                 push_value(row, row, diagonal);
			if (x < nx - 1)  push_value(row, row + 1, right);
			if (x == 0)      push_value(row, nx - 1 + nx * (y + ny * z), left);
			if (y < ny - 1)  push_value(row, row + nx, top);
		 	if (z < nz - 1)  push_value(row, row + ny * ny, front);
			if (z == 0)      push_value(row, x + ny * (y + ny * (nz - 1)), back);
		}
	}

	types::vector
	helmholtz_dirichlet::boundary(types::vector& t, types::vector& b) const
	{
		const size_type &nx = m_geometry.nx, &ny = m_geometry.ny,
		                &nz = m_geometry.nz, &n = m_geometry.n;
		types::vector r(nx * (ny - 1), nz);

		if (ny > 1) {
			for (index_type x = 0; x < nx; ++x) {
				for (index_type z = 0; z < nz; ++z) {
					const index_type i = x + nx * z;
					r[x + nx * (ny - 1) * z] = m_scale * n * n * b[i];
					r[x + nx * (ny - 2 + (ny - 1) * z)] = m_scale * n * n * t[i];
				}
			}
		}
		return r;
	}

}
