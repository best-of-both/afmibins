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

#ifndef LINALG_SPARSE_MATRIX_H
#define LINALG_SPARSE_MATRIX

#include <type_traits>

#include "types/typedefs.h"
#include "futures/futures.h"
#include "interface.h"
#include "vector.h"
#include "matrix.h"

namespace linalg {

	using types::size_type;
	using types::index_type;

	std::ostream& operator<<(std::ostream&, sparse_matrix&);

	class sparse_matrix : public base_matrix {
		public:
			typedef sparse_matrix result_type;
		private:
			size_type m_length;
		protected:
			using base_matrix::m_values;
			index_type* m_col;
			index_type* m_row_end;
			void push_value(index_type, index_type, double);
		public:
			template<typename Wrapper>
			static sparse_matrix make_from(Wrapper& w) { return sparse_matrix(w.rows(), w.cols()); }

			sparse_matrix(size_type, size_type, index_type);
			sparse_matrix(size_type, index_type);
			~sparse_matrix();

		friend void gemv(vector&, double, sparse_matrix, vector, double);
		friend std::ostream& operator<<(std::ostream&, sparse_matrix&);
	};

}

namespace futures {

	using linalg::sparse_matrix;
	using linalg::vector;

	template<>
	class mult_op<sparse_matrix, vector> {
		public:
			typedef vector result_type;
			static const char repr = '*';
		private:
			static void apply(result_type&, sparse_matrix, vector);
			static future<sparse_matrix, vector, mult_op> wrap(sparse_matrix, vector);
			static size_type cols(sparse_matrix, vector x) { return x.cols(); }
			static size_type rows(sparse_matrix a, vector) { return a.rows(); }

		template<typename, typename, typename> friend class mult_op;
		friend future<sparse_matrix, vector, mult_op> operator*<>(sparse_matrix, vector);
		friend class future<sparse_matrix, vector, mult_op>;
	};

	template<typename Left>
	class mult_op<Left, sparse_matrix, typename std::enable_if<std::is_arithmetic<Left>::value>::type> {
		public:
			typedef sparse_matrix result_type;
			static void apply(result_type&, Left, sparse_matrix);
			static const char repr = '*';
		private:
			static future<double, sparse_matrix, mult_op> wrap(Left, sparse_matrix);
			static size_type cols(Left, sparse_matrix x) { return x.cols(); }
			static size_type rows(Left, sparse_matrix x) { return x.rows(); }

		template<typename, typename, typename> friend class mult_op;
		friend future<Left, sparse_matrix, mult_op> operator*<Left, sparse_matrix>(Left, sparse_matrix);
		friend class future<Left, sparse_matrix, mult_op>;
	};

	template<>
	class add_op<sparse_matrix, sparse_matrix> {
		public:
			typedef sparse_matrix result_type;
			static void apply(result_type&, sparse_matrix, sparse_matrix); // make private ?
			static const char repr = '+';
		private:
			static future<sparse_matrix, sparse_matrix, add_op> wrap(sparse_matrix, sparse_matrix);
			static size_type cols(sparse_matrix, sparse_matrix right) { return right.cols(); }
			static size_type rows(sparse_matrix, sparse_matrix right) { return right.rows(); }

		template<typename, typename, typename> friend class add_op;
		friend future<sparse_matrix, sparse_matrix, add_op> operator+<>(sparse_matrix, sparse_matrix);
		friend class future<sparse_matrix, sparse_matrix, add_op>;
	};

	template<>
	class sub_op<sparse_matrix, sparse_matrix> {
		public:
			typedef sparse_matrix result_type;
			static void apply(result_type&, sparse_matrix, sparse_matrix); // make private ?
			static const char repr = '-';
		private:
			static future<sparse_matrix, sparse_matrix, sub_op> wrap(sparse_matrix, sparse_matrix);
			static size_type cols(sparse_matrix, sparse_matrix right) { return right.cols(); }
			static size_type rows(sparse_matrix, sparse_matrix right) { return right.rows(); }

		template<typename, typename, typename> friend class sub_op;
		friend future<sparse_matrix, sparse_matrix, sub_op> operator-<>(sparse_matrix, sparse_matrix);
		friend class future<sparse_matrix, sparse_matrix, sub_op>;
	};

	template<typename Left>
	future<double, sparse_matrix, mult_op>
	mult_op<Left, sparse_matrix, typename std::enable_if<std::is_arithmetic<Left>::value>::type>::wrap(Left alpha, sparse_matrix x)
	{
		return future<double, sparse_matrix, mult_op>(alpha, x);
	}

	template<typename Left>
	void
	mult_op<Left, sparse_matrix, typename std::enable_if<std::is_arithmetic<Left>::value>::type>::apply(result_type& y, Left alpha, sparse_matrix x)
	{
		for (std::size_t i = 0; i < x.rows() * x.cols(); ++i)
			y.m_values[i] = alpha * x.m_values[i];
	}
}

#endif
