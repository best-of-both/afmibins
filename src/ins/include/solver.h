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

#ifndef INS_SOLVER_H
#define INS_SOLVER_H

#include "types/typedefs.h"
#include "types/vector.h"
#include "projection.h"
#include "helmholtz.h"

namespace ins {

	using types::size_type;
	using types::vector;

	class solver {
		private:
			const types::geometry& m_geometry;
			const double m_tol, m_mu, m_rho, m_dt;
			projection m_projector;
			const helmholtz_periodic m_phhm, m_phhs;
			const helmholtz_dirichlet m_dhhm, m_dhhs;
			vector &m_u_curr, m_u_old, m_u_adv;
			vector &m_v_curr, m_v_old, m_v_adv;
			bool m_have_half_stepped;

			void advection(vector&, vector&);
		public:
			const vector &u, &v;

			void step(vector&, vector&, vector&, vector&);
			const vector p() { return m_phhs * m_projector.phidt; }

			solver(const types::geometry&, double, double, double, double,
					vector&, vector&);
	};

}

#endif
