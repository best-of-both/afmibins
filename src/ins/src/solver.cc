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
			vector& u0, vector& v0, vector& w0) :
		m_geometry(geometry),
		m_tol(tol), m_mu(mu), m_rho(rho), m_dt(dt), m_projector(geometry, tol),
		m_phhm(geometry, mu * dt / (2 * rho)), m_phhs(geometry, -mu * dt / (2 * rho)),
		m_dhhm(geometry, mu * dt / (2 * rho)), m_dhhs(geometry, -mu * dt / (2 * rho)),
		m_u_curr(u0), m_u_old(u0), m_u_adv(geometry.nx * geometry.ny * geometry.nz),
		m_v_curr(v0), m_v_old(v0), m_v_adv(geometry.nx * (geometry.ny - 1) * geometry.nz),
		m_w_curr(w0), m_w_old(w0), m_w_adv(geometry.nx * geometry.ny * geometry.nz),
		u(u0), v(v0), w(w0)
	{
		m_projector.project(u0, v0, w0);
	}

	void
	solver::step(vector& ut, vector& ub, vector& wt, vector& wb,
			vector& fx, vector& fy, vector& fz)
	{
		advection(m_u_curr, m_v_curr, m_w_curr);

		m_projector.boundary(ut, ub, wt, wb);
		vector bdy_u = -m_phhs.boundary(ut, ub);
		//vector bdy_v = -m_dhhs.boundary(vt, vb);
		vector bdy_w = -m_phhs.boundary(wt, wb);

		if (! m_have_half_stepped) {
			const double timestep = m_dt / 2;

			m_u_old = m_u_curr;
			m_v_old = m_v_curr;
			m_w_old = m_w_curr;

			vector rx = m_u_old +   bdy_u   + timestep * (1 / m_rho * fx - m_u_adv);
			vector ry = m_v_old + /*bdy_v*/ + timestep * (1 / m_rho * fy - m_v_adv);
			vector rz = m_w_old +   bdy_w   + timestep * (1 / m_rho * fz - m_v_adv);

			algo::conjugate_gradient(m_u_curr, m_phhs, rx, m_tol);
			algo::conjugate_gradient(m_v_curr, m_dhhs, ry, m_tol);
			algo::conjugate_gradient(m_w_curr, m_phhs, rz, m_tol);
		} else {
			const double timestep = m_dt;

			vector rx = m_phhm * m_u_old +   bdy_u   + timestep * (1 / m_rho * fx - m_u_adv);
			vector ry = m_dhhm * m_v_old + /*bdy_v*/ + timestep * (1 / m_rho * fy - m_v_adv);
			vector rz = m_phhm * m_w_old +   bdy_w   + timestep * (1 / m_rho * fz - m_v_adv);

			algo::conjugate_gradient(m_u_curr, m_phhs, rx, m_tol);
			algo::conjugate_gradient(m_v_curr, m_dhhs, ry, m_tol);
			algo::conjugate_gradient(m_w_curr, m_phhs, rz, m_tol);
		}

		m_projector.project(m_u_curr, m_v_curr, m_w_curr);
		m_have_half_stepped = ! m_have_half_stepped;
	}

	void
	solver::advection(vector& u, vector& v, vector& w)
	{
		const size_type &nx = m_geometry.nx, &ny = m_geometry.ny, &nz = m_geometry.nz;
		vector uu(nx * ny * nz), ww(nx * ny * nz), uw(nx * ny * nz);
		vector uv(nx * (ny - 1) * nz), vw(nx * (ny - 1) * nz);
		vector uvt(nx * nz), uvb(nx * nz),
		       vvt(nx * nz), vvb(nx * nz),
		       vwt(nx * nz), vwb(nx * nz);
		double q, r;

		if (ny > 1) {
			vector vv(nx * (ny - 2) * nz);
			vector vvt(nx * nz), vvb(nx * nz);

			for (index_type x = 0; x < nx; ++x) {
				for (index_type y = 1; y < ny-1; ++y) {
					for (index_type z = 0; z < nz; ++z) {
						index_type i = x + nx * (y + ny * z);
						index_type px = (x + nx - 1) % nx - x,
						           pz = (z + nz - 1) % nz - z;

						q = 0.5 * (u[i] + u[i + nx]);
						uu[i] = q * q;
						q = 0.5 * (v[i - nx * z] + v[i - nx * (z + 1)]);
						vv[i - nx * (2 * z + 1)] = q * q;
						q = 0.5 * (w[i] + w[i + nx * ny * pz]);
						ww[i] = q * q;

						q = 0.5 * (u[i + nx] + u[i]);
						r = 0.5 * (v[i - nx * z + px] + v[i - nx * z]);
						uv[i - nx * z] = q * r;

						q = 0.5 * (v[i - nx * z + nx * (ny - 1) * pz] + v[i - nx * z]);
						r = 0.5 * (w[i + nx] + w[i]);
						vw[i - nx * z] = q * r;

						q = 0.5 * (u[i + nx * nz * pz] + u[i]);
						r = 0.5 * (w[i + px] + w[i]);
						uw[i] = q * r;
					}
				}
				for (index_type z = 0; z < nz; ++z) {
					index_type ixz = x + nx * z;
					index_type px = (x + nx - 1) % nx - x,
							   pz = (z + nz - 1) % nz - z;

					q = 0.5 * v[x + nx * (ny - 1) * z];
					r = 0.5 * v[x + nx * (ny - 2 + (ny - 1) * z)];
					vvb[ixz] = q * q;
					vvt[ixz] = r * r;

					q = 0.5 * (u[x + nx * ny * z + nx] + u[x + nx * ny * z]);
					r = 0.5 * (v[x + nx * z + px] + v[x + nx * (ny - 1) * z]);
					uv[x + nx * (ny - 1) * z] = q * r;

					q = 0.5 * (v[x + nx * (ny - 1) * (z + pz)] + v[x + nx * (ny - 1) * z]);
					r = 0.5 * (w[x + nx * ny * z + nx] + w[x + nx * ny * z]);
					vw[x + nx * (ny - 1) * z] = q * r;
				}
			}
			algo::div(m_u_adv, m_geometry, uu, uv, uw, uvt, uvb);
			algo::div(m_v_adv, m_geometry, uv, vv, vw, vvt, vvb);
			algo::div(m_w_adv, m_geometry, uw, vw, ww, vwt, vwb);
		} else {
			for (index_type x = 0; x < nx; ++x) {
				for (index_type y = 0; y < ny-1; ++y) {
					for (index_type z = 0; z < nz; ++z) {
						index_type i = x + nx * ny * z;
						index_type px = (x + nx - 1) % nx - x,
						           pz = (z + nz - 1) % nz - z;
						
						q = 0.5 * (u[i] + u[i + px]);
						uu[i] = q * q;
						q = 0.5 * (w[i] + w[i + nx * ny * pz]);
						ww[i] = q * q;

						q = 0.5 * (u[i + nx] + u[i]);
						r = 0.5 * (v[i - nx * z + px] + v[i - nx * z]);
						uv[i - nx * z] = q * r;

						q = 0.5 * (v[i - nx * z + nx * (ny - 1) * pz] + v[i - nx * z]);
						r = 0.5 * (w[i + nx] + w[i]);
						vw[i - nx * z] = q * r;

						q = 0.5 * (u[i + nx * ny * pz] + u[i]);
						r = 0.5 * (w[i + px] + w[i]);
						uw[i] = q * r;
					}
				}
				for (index_type z = 0; z < nz; ++z) {
					index_type i = x + nx * ny * z;
					index_type px = (x + nx - 1) % nx - x,
					           pz = (z + nz - 1) % nz - z;

						q = 0.5 * (u[i] + u[i + px]);
						uu[i] = q * q;
						q = 0.5 * (w[i] + w[i + nx * ny * pz]);
						ww[i] = q * q;

						q = 0.5 * (u[i + nx * ny * pz] + u[i]);
						r = 0.5 * (w[i + px] + w[i]);
						uw[i] = q * r;
				}
			}

			algo::div(m_u_adv, m_geometry, uu, uv, uw, uvt, uvb);
			algo::div(m_w_adv, m_geometry, uw, vw, ww, vwt, vwb);
		}
	}
}
