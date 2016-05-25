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

#include <vector>
#include <cmath>

#include "types/typedefs.h"
#include "types/point.h"
#include "eulerian_grid.h"

namespace ib {

	types::point
	eulerian_grid::coordinates(index_type index) const
	{
		const size_type &nx = m_geometry.nx, &ny = m_geometry.ny;
		const size_type &n = m_geometry.n;
		return types::point(((double) (index % nx)) / n,
		                    ((double) (index % (nx * ny) / nx)) / n,
		                    ((double) (index / (nx * ny))) / n, index);
	}

	double
	eulerian_grid::volume(index_type index) const
	{
		const size_type &n = m_geometry.n;
		return 1. / (n * n * n);
	}

	std::vector<types::point>
	eulerian_grid::search_nearby(types::point coords) const
	{
		using std::floor;
		const size_type &nx = m_geometry.nx, &ny = m_geometry.ny,
		                &nz = m_geometry.nz, &n = m_geometry.n;

		std::vector<types::point> results;
		coords = m_geometry.fit_to_box(coords);

		for (int i = 0; i < 64; ++i) {
			const int dix = i % 4, diy = (i / 4) % 4, diz = (i / 16) % 4;
			const index_type fx = floor(coords.x * n),
			                 fy = floor(coords.y * n),
			                 fz = floor(coords.z * n);
			const index_type ix = (fx + dix + nx - 1) % nx,
			                 iy = (fy + diy),
			                 iz = (fz + diz + nz - 1) % nz;
			const double x = (fx + dix - 1.0) / n,
			             y = (fy + diy - 1.0) / n,
			             z = (fz + diz - 1.0) / n;
			types::point point(x, y, z, ix + nx * (iy - 1 + ny * iz));
			if (!m_geometry.fits_to_box(point))
				continue;
			results.push_back(point);
		}

		return results;
	}

}
