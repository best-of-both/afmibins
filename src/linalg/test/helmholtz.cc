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

class helmholtz : public sparse_matrix {
	protected:
		double m_scale;
	public:
		template<typename Wrapper>
		static helmholtz make_from(Wrapper& w) { return helmholtz(w.rows(), w.cols()); }

		helmholtz(size_type, double);
};

helmholtz::helmholtz(size_type rows, double scale) :
	sparse_matrix(rows, rows, 3), m_scale(scale)
{
	for (unsigned int i = 0; i < rows; ++i)
	{
		if (i == rows - 1) push_value(i, 0,      scale);
		if (i > 0)         push_value(i, i-1,    scale);
		                   push_value(i, i,      1 - 2.0 * scale);
		if (i < rows - 1)  push_value(i, i+1,    scale);
		if (i == 0)        push_value(i, rows-1, scale);
	}
}

int
main(void)
{
	const unsigned int N = 100;
	helmholtz L(N, 0.5);
	vector x(N, 1);
	vector y = L * x;
	for (unsigned int i = 0; i < N; ++i)
		assert(y[i] == 1);
	return 0;
}
