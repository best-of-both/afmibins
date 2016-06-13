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

#include <iostream>

#include "futures/futures.h"
#include "optimizations.h"
#include "interface.h"
#include "vector.h"
#include "dense_matrix.h"

namespace linalg {

	vector&
	operator+=(vector& v, future<future<double, dense_matrix, mult_op>, vector, mult_op> w)
	{
		auto& l = w.left;
		gemv(v, l.left, l.right, w.right, 1);
		return v;
	}

	vector&
	operator+=(vector& v, future<double, future<dense_matrix, vector, mult_op>, mult_op> w)
	{
		auto& r = w.right;
		gemv(v, w.left, r.left, r.right, 1);
		return v;
	}

	vector&
	operator+=(vector& v, future<dense_matrix, future<double, vector, mult_op>, mult_op> w)
	{
		auto& r = w.right;
		gemv(v, r.left, w.left, r.right, 1);
		return v;
	}

	vector&
	operator+=(vector& v, future<dense_matrix, vector, mult_op> w)
	{
		gemv(v, 1, w.left, w.right, 1);
		return v;
	}

	vector&
	operator-=(vector& v, future<future<double, dense_matrix, mult_op>, vector, mult_op> w)
	{
		auto& l = w.left;
		gemv(v, -l.left, l.right, w.right, 1);
		return v;
	}

	vector&
	operator-=(vector& v, future<double, future<dense_matrix, vector, mult_op>, mult_op> w)
	{
		auto& r = w.right;
		gemv(v, -w.left, r.left, r.right, 1);
		return v;
	}

	vector&
	operator-=(vector& v, future<dense_matrix, future<double, vector, mult_op>, mult_op> w)
	{
		auto& r = w.right;
		gemv(v, -r.left, w.left, r.right, 1);
		return v;
	}

	vector&
	operator-=(vector& v, future<dense_matrix, vector, mult_op> w)
	{
		gemv(v, -1, w.left, w.right, 1);
		return v;
	}

}
