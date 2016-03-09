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

#include "mg/solver.h"
#include "dt/geometry.h"
#include "dt/vector.h"

template<class T>
class test_solver : public mg::solver<T> {
	public:
		using mg::solver<T>::L;
		test_solver(double tol, double nu1, double nu2)
			: mg::solver<T>(tol, nu1, nu2) {};
};

int
main(void)
{
	typedef dt::geometry<1u, 8u, 8u, 8u> geometry;

	test_solver<geometry> solver(1e-15, 2, 2);
	dt::vector<4096> b(1);
	b[0] = 2;

	solver.solve(b);
	std::cout << solver.L * b << std::endl;
}
