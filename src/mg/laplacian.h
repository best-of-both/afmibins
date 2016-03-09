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

#ifndef MG_LAPLACIAN_H
#define MG_LAPLACIAN_H

#include "dt/vector.h"
#include "dt/matrix.h"
#include "alg/smooth.h"

namespace mg {

	template<class Geometry>
	class laplacian : public dt::matrix<Geometry::nx * Geometry::ny * Geometry::nz> {
		public:
			static constexpr unsigned int nx = Geometry::nx;
			static constexpr unsigned int ny = Geometry::ny;
			static constexpr unsigned int nz = Geometry::nz;
			static constexpr unsigned int n = Geometry::n;
		private:
			using dt::matrix<nx * ny * nz>::values;
			using dt::matrix<nx * ny * nz>::row_start;
			using dt::matrix<nx * ny * nz>::column;
		public:
			laplacian();
	};

	template<class Geometry>
	laplacian<Geometry>::laplacian()
		: dt::matrix<nx * ny * nz>()
	{
		for (auto i = 0u; i < nx * ny * nz; ++ i) {
			row_start[i] = values.size();
			auto x = i % nx, y = (i % (nx * ny)) / nx, z = i / (nx * ny);
			const double offdiagonal = n * n;
			double diagonal = -6.0 * offdiagonal,
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

			if (y == 0) {
				diagonal += bottom;
				bottom = 0;
			}
			if (y == ny - 1) {
				diagonal += top;
				top = 0;
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

			if (front != 0 && z == nz - 1) {
				values.push_back(front);
				column.push_back(x + nx * y);
			}

			if (back != 0 && z > 0) {
				values.push_back(back);
				column.push_back(i - nx * ny);
			}

			if (bottom != 0 && y > 0) {
				values.push_back(bottom);
				column.push_back(i - nx);
			}

			if (right != 0 && x == nx - 1) {
				values.push_back(right);
				column.push_back(nx * (y + ny * z));
			}

			if (left != 0 && x > 0) {
				values.push_back(left);
				column.push_back(i - 1);
			}

			if (diagonal != 0) {
				values.push_back(diagonal);
				column.push_back(i);
			}

			if (right != 0 && x < nx - 1) {
				values.push_back(right);
				column.push_back(i + 1);
			}

			if (left != 0 && x == 0) {
				values.push_back(left);
				column.push_back(nx - 1 + nx * (y + ny * z));
			}

			if (top != 0 && y < ny - 1) {
				values.push_back(top);
				column.push_back(i + nx);
			}

			if (front != 0 && z < nz - 1) {
				values.push_back(front);
				column.push_back(i + nx * ny);
			}

			if (back != 0 && z == 0) {
				values.push_back(back);
				column.push_back(x + nx * (y + ny * (nz - 1)));
			}
		}
		row_start[nx * ny * nz] = values.size();
	}
}

#endif
