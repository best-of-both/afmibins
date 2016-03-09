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

#ifndef NAVIERSTOKES_H
#define NAVIERSTOKES_H

#include <iomanip>
#include "dt/vector.h"
#include "dt/matrix.h"
#include "dt/geometry.h"
#include "projection.h"
#include "helmholtz.h"
#include "alg/cg.h"
#include "alg/div.h"


namespace ns {

	template<class Geometry>
	class solver {
		private:
			static constexpr unsigned n = Geometry::n;
			static constexpr unsigned nx = Geometry::nx;
			static constexpr unsigned ny = Geometry::ny;
			static constexpr unsigned nz = Geometry::nz;
			static constexpr unsigned size = nx * ny * nz;

			void advection(
					dt::vector<nx * ny * nz>& u,
					dt::vector<nx * (ny-1) * nz>& v,
					dt::vector<nx * ny * nz>& w);

			projection<Geometry> projector;
			helmholtz_periodic<Geometry>  phhm, phhs;
			helmholtz_dirichlet<Geometry> dhhm, dhhs;
			dt::vector<nx * ny * nz>	 &u_curr, u_old, u_adv;
			dt::vector<nx * (ny-1) * nz> &v_curr, v_old, v_adv;
			dt::vector<nx * ny * nz>	 &w_curr, w_old, w_adv;
			double mu, rho, dt, tol;
			bool have_half_stepped = false;

		public:
			solver(double mu, double rho, double k, double tol,
					dt::vector<nx * ny * nz>& u0,
					dt::vector<nx * (ny-1) * nz>& v0,
					dt::vector<nx * ny * nz>& w0)
				: projector(tol),
				  phhm(mu * k / (2 * rho)), phhs(-mu * k / (2 * rho)),
				  dhhm(mu * k / (2 * rho)), dhhs(-mu * k / (2 * rho)),
				  u_curr(u0), v_curr(v0), w_curr(w0),
				  mu(mu), rho(rho), dt(k), tol(tol)
			{
				projector.project(u_curr, v_curr, w_curr);
			}

			void step(
					dt::vector<nx * nz>&,
					dt::vector<nx * nz>&,
					dt::vector<nx * nz>&,
					dt::vector<nx * nz>&,
					dt::vector<nx * ny * nz>&,
					dt::vector<nx * (ny-1) * nz>&,
					dt::vector<nx * ny * nz>&);

			const dt::vector<nx * ny * nz>& get_u() { return u_curr; }
			const dt::vector<nx * (ny-1) * nz>& get_v() { return v_curr; }
			const dt::vector<nx * ny * nz>& get_w() { return w_curr; }
			const dt::vector<nx * ny * nz> get_p() { return phhs * projector.get_phidt(); }
	};


	template<class Geometry>
	void
	solver<Geometry>::step(
			dt::vector<nx * nz>& ut,
			dt::vector<nx * nz>& ub,
			dt::vector<nx * nz>& wt,
			dt::vector<nx * nz>& wb,
			dt::vector<nx * ny * nz>& fx,
			dt::vector<nx * (ny-1) * nz>& fy,
			dt::vector<nx * ny * nz>& fz)
	{
		//static dt::vector<nx * nz> dummy;

		advection(u_curr, v_curr, w_curr);

		if (! have_half_stepped) {
			u_old = u_curr;
			v_old = v_curr;
			w_old = w_curr;

			projector.boundary(ut, ub, wt, wb, /*dummy, dummy, */dt / 2);
			auto bdy_u = -1 * phhs.boundary(ut, ub);
			//auto bdy_v = -1 * dhhs.boundary(vt, vb);
			auto bdy_w = -1 * phhs.boundary(wt, wb);

			auto rx = u_old +   bdy_u   + (dt / 2) * (1 / rho * fx - u_adv);
			auto ry = v_old + /*bdy_v*/ + (dt / 2) * (1 / rho * fy - v_adv);
			auto rz = w_old +   bdy_w   + (dt / 2) * (1 / rho * fz - w_adv);

			alg::conjugate_gradient(u_curr, phhs, rx, tol);
			alg::conjugate_gradient(v_curr, dhhs, ry, tol);
			alg::conjugate_gradient(w_curr, phhs, rz, tol);
		} else {
			projector.boundary(ut, ub, wt, wb, /*dummy, dummy, */dt);
			auto bdy_u = phhm.boundary(ut, ub) - phhs.boundary(ut, ub);
			//auto bdy_v = dhhm.boundary(vt, vb) - dhhs.boundary(vt, vb);
			auto bdy_w = phhm.boundary(wt, wb) - phhs.boundary(wt, wb);

			auto rx = phhm * u_old +   bdy_u   + dt * (1 / rho * fx - u_adv);
			auto ry = dhhm * v_old + /*bdy_v*/ + dt * (1 / rho * fy - v_adv);
			auto rz = phhm * w_old +   bdy_w   + dt * (1 / rho * fz - w_adv);

			alg::conjugate_gradient(u_curr, phhs, rx, tol);
			alg::conjugate_gradient(v_curr, dhhs, ry, tol);
			alg::conjugate_gradient(w_curr, phhs, rz, tol);
		}

		projector.project(u_curr, v_curr, w_curr);

		have_half_stepped = ! have_half_stepped;
	}

	template<class Geometry>
	void
	solver<Geometry>::advection(
			dt::vector<nx * ny * nz>& u,
			dt::vector<nx * (ny-1) * nz>& v,
			dt::vector<nx * ny * nz>& w)
	{
		dt::vector<nx * ny * nz> uu, ww, uw;
		dt::vector<nx * (ny-1) * nz> uv, vw;
		dt::vector<nx * nz> uvt, uvb, vvt, vvb, vwt, vwb;
		double q, r;

		if (ny > 1) {
			dt::vector<nx * (ny-2) * nz> vv;
			dt::vector<nx * nz> vvt, vvb;

			for (auto i = 0u; i < nx * ny * nz; ++i) {
				auto x = i % nx,
				     y = (i % (nx * ny)) / nx,
				     z = i / (nx * ny),
				     ixz = x + nx * z;

				q = 0.5 * (u[i] + u[i + ((x + nx - 1) % nx - x)]);
				uu[i] = q * q;

				if (y == 0u) {
					/* v^2 at bottom */
					q = 0.5 * (v[i - nx * z]);
					vvb[ixz] = q * q;
				}
				else if (y == ny - 1) {
					/* v^2 at top */
					q = 0.5 * (v[i - nx * (z + 1)]);
					vvt[ixz] = q * q;
				}
				else {
					q = 0.5 * (v[i - nx * z] + v[i - nx * (z + 1)]);
					vv[i - nx * (2 * z + 1)] = q * q;
				}

				q = 0.5 * (w[i] + w[i + nx * ny * ((z + nz - 1) % nz - z)]);
				ww[i] = q * q;

				if (y < ny - 1) {
					q = 0.5 * (u[i + nx] + u[i]);
					r = 0.5 * (v[i - nx * z + ((x + nx - 1) % nx - x)] + v[i - nx * z]);
					uv[i - nx * z] = q * r;

					q = 0.5 * (v[i - nx * z + nx * (ny - 1) * ((z + nz - 1) % nz - z)] + v[i - nx * z]);
					r = 0.5 * (w[i + nx] + w[i]);
					vw[i - nx * z] = q * r;
				}

				q = 0.5 * (u[i + nx * ny * ((z + nz - 1) % nz - z)] + u[i]);
				r = 0.5 * (w[i + ((x + nx - 1) % nx - x)] + w[i]);
				uw[i] = q * r;
			}

			alg::div<nx, ny, nz, n>(u_adv, uu, uv, uw, uvt, uvb);
			alg::div<nx, ny-1, nz, n>(v_adv, uv, vv, vw, vvt, vvb);
			alg::div<nx, ny, nz, n>(w_adv, uw, vw, ww, vwt, vwb);
			return;
		} else {
			for (auto i = 0u; i < nx * ny * nz; ++i) {
				auto x = i % nx,
				     y = (i % (nx * ny)) / nx,
				     z = i / (nx * ny);

				q = 0.5 * (u[i] + u[i + ((x + nx - 1) % nx - x)]);
				uu[i] = q * q;

				q = 0.5 * (w[i] + w[i + nx * ny * ((z + nz - 1) % nz - z)]);
				ww[i] = q * q;

				if (y < ny - 1) {
					q = 0.5 * (u[i + nx] + u[i]);
					r = 0.5 * (v[i - nx * z + ((x + nx - 1) % nx - x)] + v[i - nx * z]);
					uv[i - nx * z] = q * r;

					q = 0.5 * (v[i - nx * z + nx * (ny - 1) * ((z + nz - 1) % nz - z)] + v[i - nx * z]);
					r = 0.5 * (w[i + nx] + w[i]);
					vw[i - nx * z] = q * r;
				}

				q = 0.5 * (u[i + nx * ny * ((z + nz - 1) % nz - z)] + u[i]);
				r = 0.5 * (w[i + ((x + nx - 1) % nx - x)] + w[i]);
				uw[i] = q * r;
			}

			alg::div<nx, ny, nz, n>(u_adv, uu, uv, uw, uvt, uvb);
			alg::div<nx, ny, nz, n>(w_adv, uw, vw, ww, vwt, vwb);
			return;
		}

	}

}

#endif
