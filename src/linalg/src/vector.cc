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

#include <cassert>
#include <cmath>

#include "types/typedefs.h"
#include "matrix.h"
#include "vector.h"

namespace linalg {

	using types::size_type;
	using types::index_type;

	vector::vector(size_type rows, double v) :
		base_matrix(rows, 1, new double[rows])
	{
		for (index_type i = 0; i < rows; ++i)
			m_values[i] = v;
	}

	vector::vector(std::initializer_list<double> list) :
		base_matrix(list.size(), 1, new double[list.size()])
	{
		index_type index = 0;
		for (double v: list)
			m_values[index++] = v;
	}


	vector&
	vector::operator+=(vector other)
	{
		assert(size() == other.size());
		for (std::size_t i = 0; i < size(); ++i)
			operator[](i) += other[i];
		return *this;
	}

	vector&
	vector::operator-=(vector other)
	{
		assert(size() == other.size());
		for (std::size_t i = 0; i < size(); ++i)
			operator[](i) -= other[i];
		return *this;
	}

	vector&
	vector::operator*=(double v)
	{
		for (std::size_t i = 0; i < size(); ++i)
			operator[](i) *= v;
		return *this;
	}

	double
	abs(vector v) {
		return sqrt(v * v);
	}

	std::ostream&
	operator<<(std::ostream& out, vector& v)
	{
		out << "[";
		for (std::size_t i = 0; i < v.rows(); ++i) {
			out << v[i];
			if (i < v.rows() - 1)
				out << " ";
		}
		out << "]'";
		return out;
	}

}

namespace futures {

	using linalg::vector;

	future<double, vector, mult_op>
	mult_op<double, vector>::wrap(double alpha, vector x)
	{
		return future<double, vector, mult_op>(alpha, x);
	}

	void
	mult_op<double, vector>::apply(result_type& y, double alpha, vector x)
	{
		for (std::size_t i = 0; i < x.rows(); ++i)
			y[i] = alpha * x[i];
	}

	future<vector, vector, mult_op>
	mult_op<vector, vector>::wrap(vector x, vector y)
	{
		assert(x.size() == y.size());
		return future<vector, vector, mult_op>(x, y);
	}

	void
	mult_op<vector, vector>::apply(result_type& alpha, vector x, vector y)
	{
		alpha = 0;
		for (std::size_t i = 0; i < x.size(); ++i)
			alpha += x[i] * y[i];
	}

	future<vector, vector, add_op>
	add_op<vector, vector>::wrap(vector left, vector right)
	{
		assert(left.rows() == right.rows());
		return future<vector, vector, add_op>(left, right);
	}

	void
	add_op<vector, vector>::apply(result_type& output, vector left, vector right)
	{
		for (std::size_t i = 0; i < right.rows(); ++i)
			output[i] = left[i] + right[i];
	}

	future<vector, vector, sub_op>
	sub_op<vector, vector>::wrap(vector left, vector right)
	{
		assert(left.rows() == right.rows());
		return future<vector, vector, sub_op>(left, right);
	}

	void
	sub_op<vector, vector>::apply(result_type& output, vector left, vector right)
	{
		for (std::size_t i = 0; i < right.rows(); ++i)
			output[i] = left[i] - right[i];
	}

}
