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

#ifndef INS_PROJECTION_H
#define INS_PROJECTION_H

#include "types/typedefs.h"
#include "types/vector.h"
#include "multigrid/solver.h"

namespace ins {

	using types::size_type;
	using types::vector;

	class projection {
		protected:
			const types::geometry& m_geometry;
			vector m_phidtx, m_phidty, m_phidtz, m_phidt;
			const multigrid::solver m_solver;
		public:
			const vector& phidt;

			void project(vector&, vector&, vector&);
			void boundary(vector&, vector&, vector&, vector&) const;

			projection(const types::geometry& geometry, double tol) :
				m_geometry(geometry),
				m_phidtx(geometry.nx * geometry.ny * geometry.nz),
				m_phidty(geometry.nx * (geometry.ny - 1) * geometry.nz),
				m_phidtz(geometry.nx * geometry.ny * geometry.nz),
				m_phidt(geometry.nx * geometry.ny * geometry.nz),
				m_solver(tol, 2, 2, geometry),
				phidt(m_phidt) {}
	};

}

#endif
