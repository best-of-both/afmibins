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
#include "smooth.h"

namespace algo {

	using types::vector;
	using types::sparse_square_matrix;

	void
	smooth(vector& x, const sparse_square_matrix& a, const vector& b)
	{
		std::size_t row_start = 0;
		for (std::size_t row = 0; row < a.rows(); ++row) {
			double value = 0.0;
			double diagonal = 0.0;
			for (std::size_t i = row_start; i < a.m_row_end[row]; ++i) {
				if (a.m_col[i] == row) {
					diagonal = a.m_values[i];
					continue;
				}
				value += a.m_values[i] * x[a.m_col[i]];
			}
			row_start = a.m_row_end[row];
			x[row] = (b[row] - value) / diagonal;
		}
	}

}
