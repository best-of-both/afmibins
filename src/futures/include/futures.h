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

#ifndef FUTURES_FUTURES_H
#define FUTURES_FUTURES_H

#include <cassert>
#include <ostream>
#include <type_traits>
extern "C" {
#include <cblas.h>
}
#include "types/typedefs.h"

namespace futures {

	using types::size_type;

	template<typename, typename, typename = void> class mult_op;
	template<typename, typename, typename = void> class add_op;
	template<typename, typename, typename = void> class sub_op;
	template<typename, typename, template<typename, typename, typename> class> class future;
	template<typename Left, typename Right> future<Left, Right, mult_op> operator*(Left, Right);
	template<typename Left, typename Right> future<Left, Right, add_op> operator+(Left, Right);
	template<typename Left, typename Right> future<Left, Right, sub_op> operator-(Left, Right);
	template<typename Right> future<double, Right, mult_op> operator-(Right);

	template<typename Result, typename Arg>
	class maker {
		private:
			maker() {}
			static Result make(Arg& arg) { return Result::make_from(arg); }
		template<typename, typename, template<typename, typename, typename> class> friend class future;
	};

	template<typename Left, typename Right, template<typename, typename, typename> class Op>
	class future {
		public:
			typedef typename Op<Left, Right, void>::result_type result_type;
		private:
			void resolve(result_type&);
			future(Left left, Right right) :
				left(left), right(right) {}
		public:
			Left left;
			Right right;

			size_type cols() { return Op<Left, Right, void>::cols(left, right); }
			size_type rows() { return Op<Left, Right, void>::rows(left, right); }
			operator result_type();

		friend std::ostream& operator<<(std::ostream& out, future w)
		{
			out << "(" << w.left << " " << Op<Left, Right, void>::repr << " " << w.right << ")";
			return out;	
		}

		template<typename, typename, typename> friend class mult_op;
		template<typename, typename, typename> friend class add_op;
		template<typename, typename, typename> friend class sub_op;
	};

	template<typename Left, typename Right, template<typename, typename, typename> class Op>
	void
	future<Left, Right, Op>::resolve(result_type& output)
	{
		Op<Left, Right, void>::apply(output, left, right);
	}

	template<typename Left, typename Right, template<typename, typename, typename> class Op>
	future<Left, Right, Op>::operator result_type()
	{
		result_type tmp = maker<result_type, future>::make(*this);
		resolve(tmp);
		return tmp;
	}

	template<typename Left, typename Right, typename Enable>
	class mult_op {
		public:
			typedef typename mult_op<typename Left::result_type, typename Right::result_type, void>::result_type result_type;
			static const char repr = '*';
		private:
			static void apply(result_type&, Left, Right);
			static future<Left, Right, mult_op> wrap(Left, Right);
			static size_type cols(Left, Right right) { return right.cols(); }
			static size_type rows(Left left, Right) { return left.rows(); }

		template<typename, typename, typename> friend class mult_op;
		friend future<Left, Right, mult_op> operator*<>(Left, Right);
		friend class future<Left, Right, mult_op>;
	};

	template<typename Left, typename Right, typename Enable>
	future<Left, Right, mult_op>
	mult_op<Left, Right, Enable>::wrap(Left left, Right right)
	{
		assert(left.cols() == right.rows());
		return future<Left, Right, mult_op>(left, right);
	}

	template<typename Left, typename Right, typename Enable>
	void
	mult_op<Left, Right, Enable>::apply(result_type& output, Left left, Right right)
	{
		mult_op<typename Left::result_type, typename Right::result_type, void>::apply(output, left, right);
	}

	template<typename Left, typename Right>
	class mult_op<Left, Right, typename std::enable_if<std::is_arithmetic<Left>::value>::type> {
		public:
			typedef typename Right::result_type result_type;
			static void apply(result_type&, Left, Right);
			static const char repr = '*';
		private:
			static future<Left, Right, mult_op> wrap(Left, Right);
			static size_type cols(Left, Right x) { return x.cols(); }
			static size_type rows(Left, Right x) { return x.rows(); }

		template<typename, typename, typename> friend class mult_op;
		friend future<Left, Right, mult_op> operator*<>(Left, Right);
		friend future<Left, Right, mult_op> operator-<>(Right);
		friend class future<Left, Right, mult_op>;
	};

	template<typename Left, typename Right>
	future<Left, Right, mult_op>
	mult_op<Left, Right, typename std::enable_if<std::is_arithmetic<Left>::value>::type>::wrap(Left alpha, Right x)
	{
		return future<Left, Right, mult_op>(alpha, x);
	}

	template<typename Left, typename Right>
	void
	mult_op<Left, Right, typename std::enable_if<std::is_arithmetic<Left>::value>::type>::apply(result_type& y, Left alpha, Right x)
	{
		mult_op<Left, typename Right::result_type>::apply(y, alpha, x);
	}

	template<typename Left, typename Right, typename Enable>
	class add_op {
		public:
			typedef typename add_op<typename Left::result_type, typename Right::result_type>::result_type result_type;
			static void apply(result_type&, Left, Right); // make private ?
			static const char repr = '+';
		private:
			static future<Left, Right, add_op> wrap(Left, Right);
			static size_type cols(Left, Right right) { return right.cols(); }
			static size_type rows(Left left, Right) { return left.rows(); }

		template<typename, typename, typename> friend class add_op;
		friend future<Left, Right, add_op> operator+<>(Left, Right);
		friend class future<Left, Right, add_op>;
	};

	template<typename Left, typename Right, typename Enable>
	future<Left, Right, add_op>
	add_op<Left, Right, Enable>::wrap(Left left, Right right)
	{
		assert(left.rows() == right.rows());
		return future<Left, Right, add_op>(left, right);
	}

	template<typename Left, typename Right, typename Enable>
	void
	add_op<Left, Right, Enable>::apply(result_type& output, Left left, Right right)
	{
		add_op<typename Left::result_type, typename Right::result_type>::apply(output, left, right);
	}

	template<typename Left, typename Right, typename Enable>
	class sub_op {
		public:
			typedef typename sub_op<typename Left::result_type, typename Right::result_type>::result_type result_type;
			static void apply(result_type&, Left, Right); // make private ?
			static const char repr = '-';
		private:
			static future<Left, Right, sub_op> wrap(Left, Right);
			static size_type cols(Left, Right right) { return right.cols(); }
			static size_type rows(Left left, Right) { return left.rows(); }

		template<typename, typename, typename> friend class sub_op;
		friend future<Left, Right, sub_op> operator-<>(Left, Right);
		friend class future<Left, Right, sub_op>;
	};

	template<typename Left, typename Right, typename Enable>
	future<Left, Right, sub_op>
	sub_op<Left, Right, Enable>::wrap(Left left, Right right)
	{
		assert(left.rows() == right.rows());
		return future<Left, Right, sub_op>(left, right);
	}

	template<typename Left, typename Right, typename Enable>
	void
	sub_op<Left, Right, Enable>::apply(result_type& output, Left left, Right right)
	{
		sub_op<typename Left::result_type, typename Right::result_type>::apply(output, left, right);
	}

	template<typename Left, typename Right>
	future<Left, Right, mult_op>
	operator*(Left left, Right right)
	{
		return mult_op<Left, Right, void>::wrap(left, right);
	}

	template<typename Left>
	future<double, Left, mult_op>
	operator*(Left left, double right)
	{
		return right * left;
	}

	template<typename Left>
	future<double, Left, mult_op>
	operator/(Left left, double right)
	{
		return (1 / right) * left;
	}

	template<typename Left, typename Right>
	future<Left, Right, add_op>
	operator+(Left left, Right right)
	{
		return add_op<Left, Right>::wrap(left, right);
	}

	template<typename Left, typename Right>
	future<Left, Right, sub_op>
	operator-(Left left, Right right)
	{
		return sub_op<Left, Right>::wrap(left, right);
	}

	template<typename Right>
	future<double, Right, mult_op>
	operator-(Right right)
	{
		return mult_op<double, Right>::wrap(-1.0, right);
	}

	template<typename Right>
	Right
	operator+(Right right)
	{
		return right;
	}

}

using futures::operator*;
using futures::operator/;
using futures::operator+;
using futures::operator-;

#endif
