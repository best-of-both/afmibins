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

#ifndef LINALG_DENSE_MATRIX_H
#define LINALG_DENSE_MATRIX_H

// TODO:
// * dense_matrix x dense_matrix -> dense_matrix
// * vector x dense_matrix -> vector

#include "types/typedefs.h"
#include "futures/futures.h"
#include "interface.h"
#include "vector.h"
#include "matrix.h"

namespace linalg {

	using types::size_type;
	using types::index_type;

	class dense_matrix : public base_matrix {
		public:
			typedef dense_matrix result_type;
		protected:
			index_type* m_row_starts;
			index_type* m_col;
			using base_matrix::m_values;
		public:
			template<typename Wrapper>
			static dense_matrix make_from(Wrapper&);

			dense_matrix(size_type, size_type, double);
			dense_matrix(size_type, size_type);

		friend void gemv(vector&, double, dense_matrix, vector, double);
	};

	template<typename Wrapper>
	dense_matrix
	make_from(Wrapper& w)
	{
		return dense_matrix(w.rows(), w.cols());
	}

	std::ostream& operator<<(std::ostream&, dense_matrix&);

} // namespace linalg

namespace futures {

	using linalg::dense_matrix;
	using linalg::vector;

	template<>
	class mult_op<dense_matrix, vector> {
		public:
			typedef vector result_type;
			static const char repr = '*';
		private:
			static void apply(result_type&, dense_matrix, vector); // make private ?
			static future<dense_matrix, vector, mult_op> wrap(dense_matrix, vector);
			static size_type cols(dense_matrix, vector x) { return x.cols(); }
			static size_type rows(dense_matrix a, vector) { return a.rows(); }

		template<typename, typename> friend class mult_op;
		friend future<dense_matrix, vector, mult_op> operator*<>(dense_matrix, vector);
		friend class future<dense_matrix, vector, mult_op>;
	};

	template<>
	class mult_op<double, dense_matrix> {
		public:
			typedef dense_matrix result_type;
			static void apply(result_type&, double, dense_matrix);
			static const char repr = '*';
		private:
			static future<double, dense_matrix, mult_op> wrap(double, dense_matrix);
			static size_type cols(double, dense_matrix x) { return x.cols(); }
			static size_type rows(double, dense_matrix x) { return x.rows(); }

		template<typename, typename> friend class mult_op;
		friend future<double, dense_matrix, mult_op> operator*<double, dense_matrix>(double, dense_matrix);
		friend class future<double, dense_matrix, mult_op>;
	};

	template<>
	class add_op<dense_matrix, dense_matrix> {
		public:
			typedef dense_matrix result_type;
			static void apply(result_type&, dense_matrix, dense_matrix); // make private ?
			static const char repr = '+';
		private:
			static future<dense_matrix, dense_matrix, add_op> wrap(dense_matrix, dense_matrix);
			static size_type cols(dense_matrix, dense_matrix right) { return right.cols(); }
			static size_type rows(dense_matrix, dense_matrix right) { return right.rows(); }

		template<typename, typename> friend class add_op;
		friend future<dense_matrix, dense_matrix, add_op> operator+<>(dense_matrix, dense_matrix);
		friend class future<dense_matrix, dense_matrix, add_op>;
	};

	template<>
	class sub_op<dense_matrix, dense_matrix> {
		public:
			typedef dense_matrix result_type;
			static void apply(result_type&, dense_matrix, dense_matrix); // make private ?
			static const char repr = '-';
		private:
			static future<dense_matrix, dense_matrix, sub_op> wrap(dense_matrix, dense_matrix);
			static size_type cols(dense_matrix, dense_matrix right) { return right.cols(); }
			static size_type rows(dense_matrix, dense_matrix right) { return right.rows(); }

		template<typename, typename> friend class sub_op;
		friend future<dense_matrix, dense_matrix, sub_op> operator-<>(dense_matrix, dense_matrix);
		friend class future<dense_matrix, dense_matrix, sub_op>;
	};

}

#endif
