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

#ifndef NS_PROJECTION_H
#define NS_PROJECTION_H

#include "dt/vector.h"
#include "mg/solver.h"
#include "alg/div.h"
#include "alg/grad.h"

namespace ns {

	template<class Geometry>
	class projection {
		private:
			static constexpr unsigned n = Geometry::n;
			static constexpr unsigned nx = Geometry::nx;
			static constexpr unsigned ny = Geometry::ny;
			static constexpr unsigned nz = Geometry::nz;

			dt::vector<nx * ny * nz> phidtx, phidtz, phidt;
			dt::vector<nx * (ny-1) * nz> phidty; 
			mg::solver<Geometry> solver;
		public:
			projection(double tol): solver(tol, 2u, 2u) {};
			void project(
					dt::vector<nx * ny * nz>&,
					dt::vector<nx * (ny-1) * nz>&,
					dt::vector<nx * ny * nz>&);
			void boundary(
					dt::vector<nx * nz>&,
					dt::vector<nx * nz>&,
					dt::vector<nx * nz>&,
					dt::vector<nx * nz>&,
//					dt::vector<nx * nz>&,
//					dt::vector<nx * nz>&,
					double);
			dt::vector<nx * ny * nz>& get_phidt() { return phidt; }
	};

	template<class Geometry>
	void projection<Geometry>::project(
			dt::vector<nx * ny * nz>& u,
			dt::vector<nx * (ny-1) * nz>& v,
			dt::vector<nx * ny * nz>& w)
	{
		static dt::vector<nx * nz> dummy;
		alg::div<nx, ny, nz, n>(phidt, u, v, w, dummy, dummy);
		solver.solve(phidt);
		alg::grad<nx, ny, nz, n>(phidtx, phidty, phidtz, phidt);
		u -= phidtx;
		v -= phidty;
		w -= phidtz;
	}

	template<class Geometry>
	void
	projection<Geometry>::boundary(
			dt::vector<nx * nz>& ut,
			dt::vector<nx * nz>& ub,
			dt::vector<nx * nz>& wt,
			dt::vector<nx * nz>& wb,
//			dt::vector<nx * nz>& fyt,
//			dt::vector<nx * nz>& fyb,
			double /*dt*/)
	{
		if (ny > 1) {
			for (auto i = 0u; i < nx * ny; ++i) {
				auto x = i % nx, z = i / nx,
				     ixzb = x + nx * ny * z,
				     ixzt = x + nx * (ny * (z + 1) - 1);
				double pxb1 = phidtx[ixzb], pxb2 = phidtx[ixzb + nx],
				       pxt1 = phidtx[ixzt], pxt2 = phidtx[ixzt - nx],
				       pzb1 = phidtz[ixzb], pzb2 = phidtx[ixzb + nx],
				       pzt1 = phidtz[ixzt], pzt2 = phidtx[ixzt - nx];
//				       fybxh = (fyb[i + ((x + 1) % nx - x)] - fyb[i]),
//				       fytxh = (fyt[i + ((x + 1) % nx - x)] - fyt[i]),
//				       fybzh = (fyb[i + nx * ((z + 1) % nz - z)] - fyb[i]),
//				       fytzh = (fyt[i + nx * ((z + 1) % nz - z)] - fyt[i]);

				ut[i] += (9 * pxt1 - pxt2/* - 3 * dt * fytxh*/) / 8;
				ub[i] += (9 * pxb1 - pxb2/* - 3 * dt * fybxh*/) / 8;
				wt[i] += (9 * pzt1 - pzt2/* - 3 * dt * fytzh*/) / 8;
				wb[i] += (9 * pzb1 - pzb2/* - 3 * dt * fybzh*/) / 8;
			}
		} else {
			 for (auto i = 0u; i < nx * ny; ++i) {
				auto x = i % nx, z = i / nx,
					 ixz = x + nx * ny * z;
				double px = phidtx[ixz],
				       pz = phidtz[ixz];
//				       fybxh = (fyb[i + ((x + 1) % nx - x)] - fyb[i]),
//				       fytxh = (fyt[i + ((x + 1) % nx - x)] - fyt[i]),
//				       fybzh = (fyb[i + nx * ny * ((z + 1) % nz - z)] - fyb[i]),
//				       fytzh = (fyt[i + nx * ny * ((z + 1) % nz - z)] - fyt[i]);

				ut[i] += px; // + (fybxh + fytxh) / 4;
				ub[i] += px; // - (fybxh + fytxh) / 4;
				wt[i] += pz; // + (fybzh + fytzh) / 4;
				wb[i] += pz; // - (fybzh + fytzh) / 4;
			}
		}
	}
}

#endif
