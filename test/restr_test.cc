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
class restriction_test_solver : mg::solver<T> {
	public:
		using mg::solver<T>::restriction;

		restriction_test_solver(double tol, unsigned int nu1, unsigned int nu2)
			: mg::solver<T>(tol, nu1, nu2) {}
};

int
main(void)
{
	typedef dt::geometry<1u, 1u, 1u, 1u> geometry;

	restriction_test_solver<geometry> solver(1e-16, 2, 2);

	dt::vector<8> fine = {1.4457349937345032636e-19, -1.4222019353457418616e-19,
	                      7.5757120710108818818e-20, 1.8887852354789130725e-21,
	                      -1.0404187659601848787e-19, 9.443926177393482027e-22,
	                      1.888785235478816776e-21, 2.4074124304840448163e-35};
	dt::vector<1> coarse;

	solver.restriction(coarse, fine);

	std::cout << coarse << std::endl;
}
