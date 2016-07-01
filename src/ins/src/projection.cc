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
	projection::project(vector& u, vector& v) {
		const size_type &nx = m_geometry.nx;
		const vector dummy(nx);
		algo::div(m_phidt, m_geometry, u, v, dummy, dummy);
		m_solver.solve(m_phidt);
		algo::grad(m_phidtx, m_phidty, m_geometry, m_phidt);
		u -= m_phidtx;
		v -= m_phidty;
	}

	void
	projection::boundary(vector& ut, vector& ub) const
	{
		const size_type &nx = m_geometry.nx, &ny = m_geometry.ny;
		if (ny > 1) {
			for (index_type i = 0; i < nx; ++i) {
				index_type x = i,
				           ixzb = x,
				           ixzt = x + nx * (ny - 1);
				ut[i] += (9 * m_phidtx[ixzt] - m_phidtx[ixzt - nx]) / 8;
				ub[i] += (9 * m_phidtx[ixzb] - m_phidtx[ixzb + nx]) / 8;
			}
		} else {
			for (index_type i = 0; i < nx; ++i) {
				index_type x = i % nx,
				           ixz = x + nx;
				double px = m_phidtx[ixz];
				ut[i] += px;
				ub[i] += px;
			}
		}
	}

}
