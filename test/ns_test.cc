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

#include "dt/geometry.h"
#include "ns/solver.h"

int
main(void)
{
	typedef dt::geometry<4u, 1u, 1u, 1u> geometry;
	typedef ns::solver<geometry> solver;

	constexpr unsigned int nx = geometry::nx;
	constexpr unsigned int ny = geometry::ny;
	constexpr unsigned int nz = geometry::nz;

	const double mu = 1., rho = 1., dt = 0.00005, tol=1e-15;
	const unsigned int steps = 1;

	/* x and z dimensions are periodic, the top and bottom have Dirichlet
	 * boundary conditions. Due to the MAC grid, the v components have
	 * nx * nz fewer components. The v boundary is hard-coded to be zero.
	 */

	// Flow field components and boundary values
	dt::vector<nx * ny * nz> u0, w0;
	dt::vector<nx * (ny-1) * nz> v0;
	dt::vector<nx * nz> ut, ub, wt, wb;

	// Force components
	dt::vector<nx * ny * nz> fx, fz;
	dt::vector<nx * (ny-1) * nz> fy;

	/* Here is an example to initialize u0:
	 * constexpr unsigned int n = geometry::n;
	 * const unsigned int H = ny / n;
	 * for (auto i = 0u; i < nx * ny * nz; ++i) {
	 *     auto ix = i % nx, iy = (i % (nx * ny)) / nx, iz = i / (nx * ny);
	 *     double y = (double) iy / n;
	 *
	 *     u0[i] = (H - y) * y;
	 * }
	 *
	 * Vectors also implement some of the std::vector constructors (or some
	 * semblance), e.g.
	 * dt::vector<nx * ny * nz> ones(1.0); // all ones
	 */


	solver s(mu, rho, dt, tol, u0, v0, w0);

	for (auto step = 0u; step < steps; ++step) {
		// ...

		// First call to step performs the BE step
		s.step(ut, ub, wt, wb, fx, fy, fz);
		// u0, v0, w0 are now at the half time-step

		// ...

		// Second call to step performs the CN step
		s.step(ut, ub, wt, wb, fx, fy, fz);
		// u0, v0, w0 are now at the half time-step

		// To get the pressure, e.g.,
		// auto p = s.get_p();

		// ...
	}

	return 0;
}
