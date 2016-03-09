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

#ifndef ALG_CG_H
#define ALG_CG_H

#include <iostream>
#include <cmath>
#include "dt/vector.h"
#include "dt/matrix.h"

namespace alg {

	template<unsigned int Size>
	unsigned int
	conjugate_gradient(dt::vector<Size>& x, dt::matrix<Size>& A, dt::vector<Size>& b, double tol)
	{
		auto copy = b;
		dt::vector<Size> r, p, z;
		r = p = b - A * x;
		auto n0 = r * r;

		auto count = 0u;
		while (sqrt(n0) > tol) {
			z = A * p;
			auto nu = n0 / (p * z);
			x += nu * p;
			z *= nu;
			r -= z;
			auto n1 = r * r;
			auto mu = n1 / n0;
			n0 = n1;
			p *= mu;
			p += r;
			++count;
		}
		return count;
	}

}
#endif
