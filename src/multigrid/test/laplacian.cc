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
#include <iostream>

#include "types/vector.h"
#include "types/geometry.h"
#include "laplacian.h"

using namespace multigrid;

int
main(void)
{
	const types::geometry g(1, 1, 1, 10);
	types::vector x(10 * 10 * 10, 3);
	laplacian l(g);

	types::vector y = l * x;

	for (auto v: y)
		assert(v == 0);

	return 0;
}
