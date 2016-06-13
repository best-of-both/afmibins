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
#include "futures/futures.h"
#include "sparse_matrix.h"

namespace linalg {

	using types::size_type;
	using types::index_type;

	sparse_matrix::sparse_matrix(size_type rows, size_type cols, index_type els_per_row) :
		base_matrix(rows, cols, new double[els_per_row * rows]), m_length(0),
		m_col(new index_type[els_per_row * rows]),
		m_row_end(new index_type[rows])
	{
		for (unsigned int i = 0; i < rows; ++i)
			m_row_end[i] = 0;
	}

	sparse_matrix::sparse_matrix(size_type rows, index_type els_per_row) :
		base_matrix(rows, rows, new double[els_per_row * rows]), m_length(0),
		m_col(new index_type[els_per_row * rows]),
		m_row_end(new index_type[rows])
	{
		for (unsigned int i = 0; i < rows; ++i)
			m_row_end[i] = 0;
	}

	sparse_matrix::~sparse_matrix()
	{
		if (is_last_ref()) {
			delete[] m_row_end;
			delete[] m_col;
		}
	}

	void
	sparse_matrix::push_value(index_type row, index_type col, double value)
	{
		if (value == 0) return;
		m_values[m_length] = value;
		m_col[m_length] = col;
		m_row_end[row] = ++m_length;
	}

	std::ostream&
	operator<<(std::ostream& out, sparse_matrix& m)
	{
		index_type offset = 0;
		for (index_type row = 0; row < m.rows(); ++row) {
			while (offset < m.m_row_end[row]) {
				out << "(" << row << ", " << m.m_col[offset] << ") = " << m.m_values[offset] << std::endl;
				++offset;
			}
		}
		return out;
	}

}

namespace futures {

	using linalg::sparse_matrix;
	using linalg::vector;

	future<sparse_matrix, vector, mult_op>
	mult_op<sparse_matrix, vector>::wrap(sparse_matrix A, vector x)
	{
		assert(A.cols() == x.rows());
		return future<sparse_matrix, vector, mult_op>(A, x);
	}

	void
	mult_op<sparse_matrix, vector>::apply(result_type& y, sparse_matrix A, vector x)
	{
		gemv(y, 1, A, x, 0);
	}

	future<sparse_matrix, sparse_matrix, add_op>
	add_op<sparse_matrix, sparse_matrix>::wrap(sparse_matrix left, sparse_matrix right)
	{
		assert(left.cols() == right.cols() && left.rows() == right.rows());
		return future<sparse_matrix, sparse_matrix, add_op>(left, right);
	}

	void
	add_op<sparse_matrix, sparse_matrix>::apply(result_type& output, sparse_matrix left, sparse_matrix right)
	{
		for (std::size_t i = 0; i < right.rows() * right.cols(); ++i)
			output.m_values[i] = left.m_values[i] + right.m_values[i];
	}

	future<sparse_matrix, sparse_matrix, sub_op>
	sub_op<sparse_matrix, sparse_matrix>::wrap(sparse_matrix left, sparse_matrix right)
	{
		assert(left.cols() == right.cols() && left.rows() == right.rows());
		return future<sparse_matrix, sparse_matrix, sub_op>(left, right);
	}

	void
	sub_op<sparse_matrix, sparse_matrix>::apply(result_type& output, sparse_matrix left, sparse_matrix right)
	{
		for (std::size_t i = 0; i < right.rows() * right.cols(); ++i)
			output.m_values[i] = left.m_values[i] - right.m_values[i];
	}

}
