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

#ifndef LINALG_MATRIX_H
#define LINALG_MATRIX_H

#include "types/typedefs.h"
#include "futures/futures.h"

namespace linalg {

	using types::size_type;
	using types::index_type;

	typedef struct matrix_size {
		const size_type rows, cols;
		matrix_size(size_type rows, size_type cols)
			: rows(rows), cols(cols) {}
		bool operator==(matrix_size);
		matrix_size operator*(matrix_size);
		matrix_size operator+(matrix_size);
		matrix_size operator-(matrix_size);
		matrix_size& operator=(const matrix_size& other) { assert(operator==(other)); }
		int operator+(int n) { return rows * cols + n; }
		int operator-(int n) { return rows * cols - n; }
		operator size_type() { return rows * cols; }
	} matrix_size;

	class base_matrix {
		private:
			index_type& m_refs;
		protected:
			matrix_size m_size;
			double *m_values;

			bool is_last_ref() { return m_refs == 1; }
		public:
			size_type rows() const { return m_size.rows; }
			size_type cols() const { return m_size.cols; }
			matrix_size size() const { return m_size; }

			base_matrix& operator=(const base_matrix&);

			base_matrix(size_type, size_type);
			base_matrix(size_type, size_type, double*);
			base_matrix(base_matrix&&);
			base_matrix(const base_matrix&);
			~base_matrix();

		template<typename, typename> friend class futures::mult_op;
		template<typename, typename> friend class futures::add_op;
		template<typename, typename> friend class futures::sub_op;
	};

}

#endif
