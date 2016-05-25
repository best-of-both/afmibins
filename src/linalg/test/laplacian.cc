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
#include "vector.h"
#include "sparse_matrix.h"

using namespace linalg;

class laplacian : public sparse_matrix {
	public:
		template<typename Wrapper>
		static laplacian make_from(Wrapper& w) { return laplacian(w.rows(), w.cols()); }

		laplacian(size_type);
};

laplacian::laplacian(size_type rows) :
	sparse_matrix(rows, rows, 3)
{
	for (unsigned int i = 0; i < rows; ++i)
	{
		if (i == rows - 1) push_value(i, 0, 1.0);
		if (i > 0)         push_value(i, i-1, 1.0);
		                   push_value(i, i, -2.0);
		if (i < rows - 1)  push_value(i, i+1, 1);
		if (i == 0)        push_value(i, rows-1, 1.0);
	}
}

int
main(void)
{
	const unsigned int N = 100;
	laplacian L(N);
	vector x(N, 1);
	vector y = L * x;
	for (unsigned int i = 0; i < N; ++i)
		assert(y[i] == 0);
	return 0;
}
