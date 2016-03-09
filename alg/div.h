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

#ifndef ALG_DIV_H
#define ALG_DIV_H

#include "dt/vector.h"

namespace alg {

	template<unsigned int nx, unsigned int ny, unsigned int nz, unsigned int n>
	void
	div(dt::vector<nx * ny * nz>& d,
	    dt::vector<nx * ny * nz>& u,
	    dt::vector<nx * (ny-1) * nz>& v,
	    dt::vector<nx * ny * nz>& w,
	    dt::vector<nx * nz>& top,
	    dt::vector<nx * nz>& bottom)
	{
		for (auto i = 0u; i < nx * ny * nz; ++i) {
			auto x = i % nx,
			     y = (i % (nx * ny)) / nx,
			     z = i / (nx * ny),
			     ixz = x + nx * z;
			double dv;

			if (ny > 1) {
				if (y == 0u)
					dv = v[i - nx * z] - bottom[ixz];
				else if (y == ny - 1)
					dv = top[ixz] - v[i - nx * (z + 1)];
				else
					dv = v[i - nx * z] - v[i - nx * (z + 1)];
			}
			else
				dv = top[ixz] - bottom[ixz];
			d[i] = n * (
				(u[i + ((x + 1) % nx - x)] - u[i]) + dv +
				(w[i + nx * ny * ((z + 1) % nz - z)] - w[i])
			);
		}
	}

}

#endif
