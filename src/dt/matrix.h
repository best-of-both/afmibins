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

#ifndef DT_MATRIX_H
#define DT_MATRIX_H

#include <vector>
#include <utility>
#include <ostream>
#include <iostream>
#include <iomanip>
#include "vector.h"
#include "alg/smooth.h"

namespace dt {

	template<unsigned int Size>
	class matrix {
		public:
			typedef double value_type;
		private:
			typedef std::vector<std::size_t> index_vector_type;
			typedef std::vector<value_type> vector_type;
		public:
			typedef std::size_t size_type;
		protected:
			index_vector_type row_start;
			index_vector_type column;
			vector_type values;
		public:
			matrix()
				: row_start(Size + 1, 0) { };
			template<template<unsigned int> class V>
				V<Size> operator*(V<Size>);
			size_type size() { return Size; }

		friend std::ostream& operator<<(std::ostream& out, const matrix<Size>& mat)
		{
			//out << "[";
			for (std::size_t row = 0; row < Size; ++row)
			{
				std::size_t end = mat.row_start[row + 1];
				std::size_t offset = mat.row_start[row];
				for (std::size_t col = 0; col < Size; ++col)
				{
					if (end > offset && col == mat.column[offset]) {
						// out << std::setw(8) << mat.values[offset];
						out << '(' << (row + 1) << ", " << (col + 1) << ") = " << mat.values[offset] << std::endl;
						offset++;
					} // else out << std::setw(8) << "0";
					// if (col != Size - 1) out << " ";
				}
				// if (row != n - 1) out << "\n ";
			}
			// out << "]";
			return out;
		}


		friend void alg::smooth<>(dt::vector<Size>&, const matrix<Size>&, const dt::vector<Size>&);
	};

	template<unsigned int Size> template<template<unsigned int> class V>
	V<Size>
	matrix<Size>::operator*(V<Size> right)
	{
		std::size_t row = 0;
		V<Size> result;
		for (std::size_t i = 0; i < values.size(); ++i)
		{
			while (row < Size - 1 && i >= row_start[row+1]) { ++row; }
			result[row] += values[i] * right[column[i]];
		}

		return result;
	}

}

#endif
