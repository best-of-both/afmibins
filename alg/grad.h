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

#ifndef ALG_GRAD_H
#define ALG_GRAD_H

#include "dt/vector.h"

namespace alg {

	template<unsigned nx, unsigned ny, unsigned nz, unsigned n>
	void
	grad(dt::vector<nx * ny * nz>& px,
	     dt::vector<nx * (ny-1) * nz>& py,
	     dt::vector<nx * ny * nz>& pz,
	     dt::vector<nx * ny * nz>& p)
	{
		for (auto i = 0u; i < nx * ny * nz; ++i) {
			auto x = i % nx, y = (i % (nx * ny)) / nx, z = i / (nx * ny);

			px[i] = n * (p[i] - p[i + ((x + nx - 1) % nx - x)]);
			pz[i] = n * (p[i] - p[i + nx * ny * ((z + nz - 1) % nz - z)]);
			if (y > 0u)
				py[i - nx * (z + 1)] = n * (p[i] - p[i - nx]);
		}
	}

}

#endif
