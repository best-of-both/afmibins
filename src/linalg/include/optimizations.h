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

	vector& operator+=(vector&, future<future<double, dense_matrix, mult_op>, vector, mult_op>);
	vector& operator+=(vector&, future<double, future<dense_matrix, vector, mult_op>, mult_op>);
	vector& operator+=(vector&, future<dense_matrix, future<double, vector, mult_op>, mult_op>);
	vector& operator+=(vector&, future<dense_matrix, vector, mult_op>);
	vector& operator-=(vector&, future<future<double, dense_matrix, mult_op>, vector, mult_op>);
	vector& operator-=(vector&, future<double, future<dense_matrix, vector, mult_op>, mult_op>);
	vector& operator-=(vector&, future<dense_matrix, future<double, vector, mult_op>, mult_op>);
	vector& operator-=(vector&, future<dense_matrix, vector, mult_op>);

}

#endif
