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

#include "types/typedefs.h"
#include "matrix.h"

namespace linalg {

	using types::size_type;

	bool
	matrix_size::operator==(matrix_size other) {
		return other.rows == rows &&
			   other.cols == cols;
	}

	matrix_size
	matrix_size::operator*(matrix_size other) {
		assert(other.rows == cols);
		return {rows, other.cols};
	}

	matrix_size
	matrix_size::operator+(matrix_size other) {
		assert(*this == other);
		return {rows, cols};
	}

	matrix_size
	matrix_size::operator-(matrix_size other) {
		assert(*this == other);
		return {rows, cols};
	}

	base_matrix::base_matrix(size_type rows, size_type cols) :
		m_refs(*new index_type(1)), m_size(rows, cols), m_values(nullptr) {}

	base_matrix::base_matrix(size_type rows, size_type cols, double* values) :
		m_refs(*new index_type(1)), m_size(rows, cols), m_values(values) {}

	base_matrix::base_matrix(const base_matrix& other) :
		m_refs(other.m_refs), m_size(other.m_size), m_values(other.m_values)
	{
		++m_refs;
	}

	base_matrix::base_matrix(base_matrix&& other) :
		m_refs(other.m_refs), m_size(other.m_size), m_values(other.m_values)
	{
		other.m_values = nullptr;
		other.m_refs = 0;
	}

	base_matrix&
	base_matrix::operator=(const base_matrix& other)
	{
		for (unsigned int i = 0; i < size(); ++i)
			m_values[i] = other.m_values[i];
		return *this;
	}

	base_matrix::~base_matrix()
	{
		if (m_values != nullptr) {
			if (is_last_ref()) {
				delete[] m_values;
				delete &m_refs;
			}
		}
	}

}
