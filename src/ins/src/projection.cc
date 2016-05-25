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
#include "types/vector.h"
#include "algo/grad.h"
#include "algo/div.h"
#include "projection.h"

namespace ins {

	using types::size_type;
	using types::index_type;

	void
	projection::project(vector& u, vector& v, vector& w) {
		const size_type &nx = m_geometry.nx, &nz = m_geometry.nz;
		const vector dummy(nx * nz);
		algo::div(m_phidt, m_geometry, u, v, w, dummy, dummy);
		m_solver.solve(m_phidt);
		algo::grad(m_phidtx, m_phidty, m_phidtz, m_geometry, m_phidt);
		u -= m_phidtx;
		v -= m_phidty;
		w -= m_phidtz;
	}

	void
	projection::boundary(vector& ut, vector& ub, vector& wt, vector& wb) const
	{
		const size_type &nx = m_geometry.nx, &ny = m_geometry.ny, &nz = m_geometry.nz;
		if (ny > 1) {
			for (index_type i = 0; i < nx * nz; ++i) {
				index_type x = i % nx, z = i / nx,
				           ixzb = x + nx * ny * z,
				           ixzt = x + nx * (ny * (z + 1) - 1);
				ut[i] += (9 * m_phidtx[ixzt] - m_phidtx[ixzt - nx]) / 8;
				ub[i] += (9 * m_phidtx[ixzb] - m_phidtx[ixzb + nx]) / 8;
				wt[i] += (9 * m_phidtz[ixzt] - m_phidtz[ixzt - nx]) / 8;
				wb[i] += (9 * m_phidtz[ixzb] - m_phidtz[ixzb + nx]) / 8;
			}
		} else {
			for (index_type i = 0; i < nx * nz; ++i) {
				index_type x = i % nx, z = i / nx,
				           ixz = x + nx * ny * z;
				double px = m_phidtx[ixz], pz = m_phidtz[ixz];
				ut[i] += px;
				ub[i] += px;
				wt[i] += pz;
				wb[i] += pz;
			}
		}
	}

}
