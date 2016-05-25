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

#include "types/vector.h"
#include "types/matrix.h"
#include "cg.h"

namespace algo {

	using types::sparse_square_matrix;

	unsigned int
	conjugate_gradient(vector& x, const sparse_square_matrix& A, const vector& b, const double tol)
	{
		vector r = b - A * x, p = r;
		double n0 = r * r;
		double nu, n1, mu;
		const double tol2 = tol * tol;

		auto count = 0u;
		while (n0 > tol2) {
			vector z = A * p;
			nu = n0 / (p * z);
			x += nu * p;
			z *= nu;
			r -= z;
			n1 = r * r;
			mu = n1 / n0;
			n0 = n1;
			p = mu * p + r;
			++count;
		}
		return count;
	}

}
