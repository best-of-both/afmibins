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

#ifndef LINALG_OPTIMIZATIONS_H
#define LINALG_OPTIMIZATIONS_H

#include "futures/futures.h"
#include "vector.h"
#include "dense_matrix.h"

namespace linalg {

	using futures::future;
	using futures::mult_op;

	template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	vector& operator+=(vector& y, future<future<T, dense_matrix, mult_op>, vector, mult_op> f)
	{
		auto& l = f.left;
		gemv(y, l.left, l.right, f.right, 1);
		return y;
	}

	template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	vector& operator+=(vector& y, future<T, future<dense_matrix, vector, mult_op>, mult_op> f)
	{
		auto& r = f.right;
		gemv(y, f.left, r.left, r.right, 1);
		return y;
	}

	template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	vector& operator+=(vector& y, future<dense_matrix, future<T, vector, mult_op>, mult_op> f)
	{
		auto& r = f.right;
		gemv(y, r.left, f.left, r.right, 1);
		return y;
	}

	template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	vector& operator+=(vector& y, future<dense_matrix, vector, mult_op> f)
	{
		gemv(y, 1, f.left, f.right, 1);
		return y;
	}

	template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	vector& operator-=(vector& y, future<future<T, dense_matrix, mult_op>, vector, mult_op> f)
	{
		auto& l = f.left;
		gemv(y, -l.left, l.right, f.right, 1);
		return y;
	}

	template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	vector& operator-=(vector& y, future<T, future<dense_matrix, vector, mult_op>, mult_op> f)
	{
		auto& r = f.right;
		gemv(y, -f.left, r.left, r.right, 1);
		return y;
	}

	template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	vector& operator-=(vector& y, future<dense_matrix, future<T, vector, mult_op>, mult_op> f)
	{
		auto& r = f.right;
		gemv(y, -r.left, f.left, r.right, 1);
		return y;
	}
	template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	vector& operator-=(vector& y, future<dense_matrix, vector, mult_op> f)
	{
		gemv(y, -1, f.left, f.right, 1);
		return y;
	}
}

#endif
