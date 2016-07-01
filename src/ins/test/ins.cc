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

#include "types/typedefs.h"
#include "types/geometry.h"
#include "types/vector.h"
#include "ib/eulerian_grid.h"
#include "solver.h"

using namespace ins;

int
main(void)
{
	const types::geometry geometry(1, 1, 16);
	const types::size_type &nx = geometry.nx, &ny = geometry.ny,
	                       &n = geometry.n;

	const double mu = 1., rho = 1., dt = 0.00001, tol=1e-15;
	const unsigned int steps = 1;

	/* x and z dimensions are periodic, the top and bottom have Dirichlet
	 * boundary conditions. Due to the MAC grid, the v components have
	 * nx * nz fewer components. The v boundary is hard-coded to be zero.
	 */

	// Flow field components and boundary values
	ib::eulerian_grid u0(geometry), v0(geometry, nx * (ny - 1));
	types::vector ut(nx), ub(nx);

	// Force components
	ib::eulerian_grid fx(geometry, 1.0), fy(geometry, nx * (ny - 1));

	/* Here is an example to initialize u0: */
	 const types::size_type &H = geometry.height;
	 for (auto i = 0u; i < nx * ny; ++i) {
	      auto iy = i / nx;
	      double y = (iy + 0.5) / n;
	
	      u0[i] = (H - y) * y;
	 }

	 /*
	 * Vectors also implement some of the std::vector constructors (or some
	 * semblance), e.g.
	 * dt::vector<nx * ny * nz> ones(1.0); // all ones
	 */


	solver s(geometry, mu, rho, dt, tol, u0, v0);
	//std::cout << u0 << std::endl;

	for (auto step = 0u; step < steps; ++step) {
		// ...

		// First call to step performs the BE step
		s.step(ut, ub, fx, fy);
		// u0, v0, w0 are now at the half time-step

		// fx.interpolate(lagrangian_grid_inst_x)
		// fy.interpolate(lagrangian_grid_inst_y)

		// Second call to step performs the CN step
		s.step(ut, ub, fx, fy);
		// u0, v0, w0 are now at the full time-step

		// fx.interpolate(lagrangian_grid_inst_x)
		// fy.interpolate(lagrangian_grid_inst_y)

		// To get the pressure, e.g.,
		// auto p = s.get_p();

		// ...
	}

	//std::cout << u0 << std::endl;


	return 0;
}
