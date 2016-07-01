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
#include "algo/div.h"
#include "algo/grad.h"
#include "algo/cg.h"
#include "projection.h"
#include "helmholtz.h"
#include "solver.h"

namespace ins {

	using types::size_type;
	using types::index_type;
	using types::vector;

	solver::solver(const types::geometry& geometry,
			double tol, double mu, double rho, double dt,
			vector& u0, vector& v0) :
		m_geometry(geometry),
		m_tol(tol), m_mu(mu), m_rho(rho), m_dt(dt), m_projector(geometry, tol),
		m_phhm(geometry, mu * dt / (2 * rho)), m_phhs(geometry, -mu * dt / (2 * rho)),
		m_dhhm(geometry, mu * dt / (2 * rho)), m_dhhs(geometry, -mu * dt / (2 * rho)),
		m_u_curr(u0), m_u_old(u0), m_u_adv(geometry.nx * geometry.ny),
		m_v_curr(v0), m_v_old(v0), m_v_adv(geometry.nx * (geometry.ny - 1)),
		u(u0), v(v0)
	{
		m_projector.project(u0, v0);
	}

	void
	solver::step(vector& ut, vector& ub, vector& fx, vector& fy)
	{
		advection(m_u_curr, m_v_curr);

		m_projector.boundary(ut, ub);
		vector bdy_u = -m_phhs.boundary(ut, ub);
		//vector bdy_v = -m_dhhs.boundary(vt, vb);

		if (! m_have_half_stepped) {
			const double timestep = m_dt / 2;

			m_u_old = m_u_curr;
			m_v_old = m_v_curr;

			vector rx = m_u_old +   bdy_u   + timestep * (fx / m_rho - m_u_adv);
			vector ry = m_v_old + /*bdy_v*/ + timestep * (fy / m_rho - m_v_adv);

			algo::conjugate_gradient(m_u_curr, m_phhs, rx, m_tol);
			algo::conjugate_gradient(m_v_curr, m_dhhs, ry, m_tol);
		} else {
			const double timestep = m_dt;

			vector rx = m_phhm * m_u_old +   bdy_u   + timestep * (fx / m_rho - m_u_adv);
			vector ry = m_dhhm * m_v_old + /*bdy_v*/ + timestep * (fy / m_rho - m_v_adv);

			algo::conjugate_gradient(m_u_curr, m_phhs, rx, m_tol);
			algo::conjugate_gradient(m_v_curr, m_dhhs, ry, m_tol);
		}

		m_projector.project(m_u_curr, m_v_curr);
		m_have_half_stepped = ! m_have_half_stepped;
	}

	void
	solver::advection(vector& u, vector& v)
	{
		const size_type &nx = m_geometry.nx, &ny = m_geometry.ny;
		vector uu(nx * ny);
		vector uv(nx * (ny - 1));
		vector uvt(nx), uvb(nx);
		double q, r;

		if (ny > 1) {
			vector vv(nx * (ny - 2));
			vector vvt(nx), vvb(nx);

			for (index_type x = 0; x < nx; ++x) {
				index_type px = (x + nx - 1) % nx - x;
				for (index_type y = 1; y < ny-1; ++y) {
					index_type i = x + nx * y;

					q = 0.5 * (u[i] + u[i + px]);
					uu[i] = q * q;
					q = 0.5 * (v[i] + v[i - nx]);
					vv[i - nx] = q * q;

					q = 0.5 * (u[i + nx] + u[i]);
					r = 0.5 * (v[i + px] + v[i]);
					uv[i] = q * r;
				}
				index_type ixz = x;

				q = 0.5 * v[x];
				r = 0.5 * v[x + nx * (ny - 2)];
				vvb[ixz] = q * q;
				vvt[ixz] = r * r;

				q = 0.5 * (u[x + nx] + u[x]);
				r = 0.5 * (v[x + px] + v[x]);
				uv[x] = q * r;
			}
			algo::div(m_u_adv, m_geometry, uu, uv, uvt, uvb);
			algo::div(m_v_adv, m_geometry, uv, vv, vvt, vvb);
		} else {
			for (index_type x = 0; x < nx; ++x) {
				index_type i = x;
				index_type px = (x + nx - 1) % nx - x;
				
				q = 0.5 * (u[i] + u[i + px]);
				uu[i] = q * q;

				q = 0.5 * (u[i + nx] + u[i]);
				r = 0.5 * (v[i + px] + v[i]);
				uv[i] = q * r;
			}
			algo::div(m_u_adv, m_geometry, uu, uv, uvt, uvb);
		}
	}
}
