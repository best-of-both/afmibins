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

#include <iostream>
#include <cmath>
#include "vector.h"
#include "sparse_matrix.h"

using namespace linalg;

double dot(vector& a, vector& b)
{
	assert(a.size() == b.size());
	double v = 0.0;
	for (unsigned int i = 0; i < a.size(); ++i)
		v += a[i] * b[i];
	return v;
}

unsigned int
conjugate_gradient(vector& x, sparse_matrix& A, vector& b, double tol)
{
	vector r = b - A * x, p = r;
	double n0 = dot(r, r);
	double nu, n1, mu;

	auto count = 0u;
	while (sqrt(n0) > tol) {
		vector z = A * p;
		nu = n0 / dot(p, z);
		x += nu * p;
		z *= nu;
		r -= z;
		n1 = dot(r, r);
		mu = n1 / n0;
		n0 = n1;
		p = mu * p + r;
		++count;
	}
	return count;
}

class identity : public sparse_matrix {
	public:
		identity(size_type rows) :
			sparse_matrix(rows, rows, 1)
	{
		for (unsigned int i = 0; i < rows; ++i)
			push_value(i, i, 1);
	}
};

int
main(void)
{
	const unsigned int N = 100;
	vector x(N, 1.01);
	identity a(N);
	vector b(N, 1);

	unsigned int iterates = conjugate_gradient(x, a, b, 0.0001);
	double v = 0.0;
	for (unsigned int i = 0; i < N; ++i)
		v += (x[i] - b[i]) * (x[i] - b[i]);
	assert(sqrt(v) < 0.0001);
	assert(iterates == 1);
	return 0;
}
