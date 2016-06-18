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

#include <mkl.h>

#include "types/typedefs.h"
#include "interface.h"
#include "vector.h"
#include "dense_matrix.h"
#include "sparse_matrix.h"

namespace linalg {

	void
	gemv(vector& y, double alpha, dense_matrix a, vector x, double beta)
	{
		cblas_dgemv(CblasRowMajor, CblasNoTrans, a.rows(), a.cols(), alpha,
				a.m_values, a.cols(), x.m_values, x.cols(), beta, y.m_values, y.cols());
	}

	void
	gemv(vector& y, double alpha, sparse_matrix a, vector x, double beta)
	{
		using types::index_type;

		index_type row = 0;
		double row_value = y[row];
		double mult_value = 0;
		for (index_type i = 0; i < a.m_length; ++i) {
			while (i >= a.m_row_end[row]) {
				y[row] = alpha * mult_value + beta * row_value;
				mult_value = 0;
				row_value = y[++row];
			}
			mult_value += a.m_values[i] * x[a.m_col[i]];
		}
		y[row] = alpha * mult_value + beta * row_value;
	}

}
