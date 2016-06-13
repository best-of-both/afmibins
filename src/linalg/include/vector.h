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

#ifndef LINALG_VECTOR_H
#define LINALG_VECTOR_H

#include <ostream>
#include <cstddef>
#include <initializer_list>
#include <type_traits>

#include "types/typedefs.h"
#include "futures/futures.h"
#include "interface.h"
#include "matrix.h"

namespace linalg {

	using types::size_type;
	using types::index_type;

	class sparse_matrix;
	class dense_matrix;

	class vector : public base_matrix {
		public:
			/* iterator typedefs */
			typedef double value_type;
			typedef std::ptrdiff_t difference_type;
			typedef types::size_type size_type;
			typedef double* pointer;
			typedef double& reference;
			typedef pointer iterator;
			typedef const double* const_pointer;
			typedef const double& const_reference;
			typedef const_pointer const_iterator;

			/* futures typedefs */
			typedef vector result_type;
		protected:
			using base_matrix::m_values;
		public:
			double& operator[](index_type n) { return m_values[n]; }
			double& operator[](index_type n) const { return m_values[n]; }

			template<typename Wrapper>
			static vector make_from(Wrapper& w) { return vector(w.rows()); }

			iterator begin() { return m_values; }
			iterator end() { return &m_values[rows() * cols()]; }
			const_iterator begin() const { return m_values; }
			const_iterator end() const { return &m_values[rows() * cols()]; }

			vector& operator+=(vector);
			vector& operator-=(vector);
			vector& operator*=(double);

			vector(size_type, double = 0.0);
			vector(std::initializer_list<double>);

		friend void gemv(vector&, double, sparse_matrix, vector, double);
		friend void gemv(vector&, double, dense_matrix, vector, double);
	};

	std::ostream& operator<<(std::ostream&, vector&);
	double abs(vector);

}

namespace futures {

	using linalg::vector;
	using types::size_type;

	template<typename Left>
	class mult_op<Left, vector, typename std::enable_if<std::is_arithmetic<Left>::value>::type> {
		public:
			typedef vector result_type;
			static void apply(result_type&, Left, vector);
			static const char repr = '*';
		private:
			static future<Left, vector, mult_op> wrap(Left, vector);
			static size_type cols(Left, vector x) { return x.cols(); }
			static size_type rows(Left, vector x) { return x.rows(); }

		template<typename, typename, typename> friend class mult_op;
		friend future<Left, vector, mult_op> operator*<Left, vector>(Left, vector);
		friend future<Left, vector, mult_op> operator-<vector>(vector);
		friend class future<Left, vector, mult_op>;
	};

	template<>
	class mult_op<vector, vector> {
		public:
			typedef double result_type;
			static void apply(result_type&, vector, vector);
			static const char repr = '*';
		private:
			static future<vector, vector, mult_op> wrap(vector, vector);
			static size_type cols(vector, vector) { return 1; }
			static size_type rows(vector, vector) { return 1; }

		template<typename, typename, typename> friend class mult_op;
		friend future<vector, vector, mult_op> operator*<vector, vector>(vector, vector);
		friend class future<vector, vector, mult_op>;
	};

	template<typename Arg>
	class maker<double, Arg> {
		private:
			maker() {}
			static double make(Arg&) { return 0.0; }
		template<typename, typename, template<typename, typename, typename> class> friend class future;
	};

	template<>
	class add_op<vector, vector> {
		public:
			typedef vector result_type;
			static void apply(result_type&, vector, vector); // make private ?
			static const char repr = '+';
		private:
			static future<vector, vector, add_op> wrap(vector, vector);
			static size_type cols(vector, vector right) { return right.cols(); }
			static size_type rows(vector, vector right) { return right.rows(); }

		template<typename, typename, typename> friend class add_op;
		friend future<vector, vector, add_op> operator+<>(vector, vector);
		friend class future<vector, vector, add_op>;
	};

	template<>
	class sub_op<vector, vector> {
		public:
			typedef vector result_type;
			static void apply(result_type&, vector, vector); // make private ?
			static const char repr = '-';
		private:
			static future<vector, vector, sub_op> wrap(vector, vector);
			static size_type cols(vector, vector right) { return right.cols(); }
			static size_type rows(vector, vector right) { return right.rows(); }

		template<typename, typename, typename> friend class sub_op;
		friend future<vector, vector, sub_op> operator-<>(vector, vector);
		friend class future<vector, vector, sub_op>;
	};

	template<typename Left>
	future<Left, vector, mult_op>
	mult_op<Left, vector, typename std::enable_if<std::is_arithmetic<Left>::value>::type>::wrap(Left alpha, vector x)
	{
		return future<Left, vector, mult_op>(alpha, x);
	}

	template<typename Left>
	void
	mult_op<Left, vector, typename std::enable_if<std::is_arithmetic<Left>::value>::type>::apply(result_type& y, Left alpha, vector x)
	{
		for (std::size_t i = 0; i < x.rows(); ++i)
			y[i] = alpha * x[i];
	}

}

#endif
