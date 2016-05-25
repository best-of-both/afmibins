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

#ifndef TYPES_MATRIX_H
#define TYPES_MATRIX_H

#include <iostream>

#include "linalg/matrix.h"
#include "linalg/sparse_matrix.h"
#include "linalg/dense_matrix.h"
#include "vector.h"

namespace types { class sparse_square_matrix; }

namespace algo {

	using types::vector;
	using types::sparse_square_matrix;

	void smooth(vector&, const sparse_square_matrix&, const vector&);
}

namespace types {

	using linalg::dense_matrix;
	using linalg::sparse_matrix;

	class sparse_square_matrix : public sparse_matrix {
		protected:
			using sparse_matrix::m_values;
			using sparse_matrix::m_col;
			using sparse_matrix::m_row_end;
		public:
			sparse_square_matrix(size_type rows, index_type els_per_row) :
				sparse_matrix(rows, rows, els_per_row) {}

		friend void algo::smooth(vector&, const sparse_square_matrix&, const vector&);
	};
}

#endif
