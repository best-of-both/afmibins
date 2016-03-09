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

#ifndef MG_MULTIGRID_H
#define MG_MULTIGRID_H

#include <vector>
#include <cmath>
#include "dt/geometry.h"
#include "dt/vector.h"
#include "laplacian.h"
#include "alg/cg.h"


namespace mg {

	template<class Geometry, unsigned Level = Geometry::level>
	class solver {
		public:
			static constexpr unsigned n = Geometry::n;
			static constexpr unsigned nx = Geometry::nx;
			static constexpr unsigned ny = Geometry::ny;
			static constexpr unsigned nz = Geometry::nz;

			typedef Geometry geometry_type;
			typedef laplacian<geometry_type> laplacian_type;
		private:
			solver<typename geometry_type::coarse> coarser;
			const double tol;
			const unsigned int nu1, nu2;
		protected:
			static laplacian_type L;
			static void interpolation(
					dt::vector<nx * ny * nz>&,
					dt::vector<nx * ny * nz / 8>&);
			static void restriction(
					dt::vector<nx * ny * nz / 8>&,
					dt::vector<nx * ny * nz>&);
			void nested_iteration(dt::vector<nx * ny * nz>&);
			void recurse(dt::vector<nx * ny * nz>&);
		public:
			unsigned int solve(dt::vector<nx * ny * nz>&);
			solver(double tol, unsigned int nu1, unsigned int nu2)
				: coarser(tol, nu1, nu2), tol(tol), nu1(nu1), nu2(nu2)  {};

			friend class solver<typename geometry_type::fine>;
			friend class solver<typename geometry_type::coarse>;
	};

	template<class Geometry, unsigned Level>
		laplacian<Geometry> solver<Geometry, Level>::L;

	template<class Geometry, unsigned Level>
	void
	solver<Geometry, Level>::interpolation(
		dt::vector<nx * ny * nz>& output,
		dt::vector<nx * ny * nz / 8>& input)
	{
		for (auto xcc = 0u; xcc < nx / 2; ++xcc) {
			for (auto ycc = 0u; ycc < ny / 2 - 1; ++ycc) {
				for (auto zcc = 0u; zcc < nz / 2; ++zcc) {
					auto xnc = (xcc + 1) % (nx / 2),
					     ync =  ycc + 1,
					     znc = (zcc + 1) % (nz / 2);
					auto xcf = (2 * xcc + 1),
					     ycf = (2 * ycc + 1),
					     zcf = (2 * zcc + 1);
					auto xnf = (xcf + 1) % nx,
					     ynf = ycf + 1,
					     znf = (zcf + 1) % nz;
					auto blbc = input[xcc + (nx / 2) * (ycc + (ny / 2) * zcc)],
					     brbc = input[xnc + (nx / 2) * (ycc + (ny / 2) * zcc)],
					     ulbc = input[xcc + (nx / 2) * (ync + (ny / 2) * zcc)],
					     urbc = input[xnc + (nx / 2) * (ync + (ny / 2) * zcc)],
					     blfc = input[xcc + (nx / 2) * (ycc + (ny / 2) * znc)],
					     brfc = input[xnc + (nx / 2) * (ycc + (ny / 2) * znc)],
					     ulfc = input[xcc + (nx / 2) * (ync + (ny / 2) * znc)],
					     urfc = input[xnc + (nx / 2) * (ync + (ny / 2) * znc)];
					auto &blbf = output[xcf + nx * (ycf + ny * zcf)],
					     &brbf = output[xnf + nx * (ycf + ny * zcf)],
					     &ulbf = output[xcf + nx * (ynf + ny * zcf)],
					     &urbf = output[xnf + nx * (ynf + ny * zcf)],
					     &blff = output[xcf + nx * (ycf + ny * znf)],
					     &brff = output[xnf + nx * (ycf + ny * znf)],
					     &ulff = output[xcf + nx * (ynf + ny * znf)],
					     &urff = output[xnf + nx * (ynf + ny * znf)];

					blbf = 27. / 64. * blbc + 9. / 64. * (brbc + ulbc + blfc) + 3. / 64. * (urbc + ulfc + brfc) + 1. / 64. * urfc;
					brbf = 27. / 64. * brbc + 9. / 64. * (urbc + blbc + brfc) + 3. / 64. * (ulbc + blfc + urfc) + 1. / 64. * ulfc;
					ulbf = 27. / 64. * ulbc + 9. / 64. * (urbc + blbc + ulfc) + 3. / 64. * (brbc + blfc + urfc) + 1. / 64. * brfc;
					urbf = 27. / 64. * urbc + 9. / 64. * (brbc + ulbc + urfc) + 3. / 64. * (blbc + ulfc + brfc) + 1. / 64. * blfc;
					blff = 27. / 64. * blfc + 9. / 64. * (brfc + ulfc + blbc) + 3. / 64. * (urfc + ulbc + brbc) + 1. / 64. * urbc;
					brff = 27. / 64. * brfc + 9. / 64. * (urfc + blfc + brbc) + 3. / 64. * (ulfc + blbc + urbc) + 1. / 64. * ulbc;
					ulff = 27. / 64. * ulfc + 9. / 64. * (urfc + blfc + ulbc) + 3. / 64. * (brfc + blbc + urbc) + 1. / 64. * brbc;
					urff = 27. / 64. * urfc + 9. / 64. * (brfc + ulfc + urbc) + 3. / 64. * (blfc + ulbc + brbc) + 1. / 64. * blbc;
				}
			}

			for (auto zcc = 0u; zcc < nz / 2; ++zcc) {
				auto xnc = (xcc + 1) % (nx / 2),
				     znc = (zcc + 1) % (nz / 2);
				auto xcf = 2 * xcc + 1,
					 zcf = 2 * zcc + 1;
				auto xnf = (xcf + 1) % nx,
				     znf = (zcf + 1) % nz;
				auto blbc = input[xcc + (nx / 2) * (ny / 2) * zcc],
					 brbc = input[xnc + (nx / 2) * (ny / 2) * zcc],
					 blfc = input[xcc + (nx / 2) * (ny / 2) * znc],
					 brfc = input[xnc + (nx / 2) * (ny / 2) * znc],
					 ulbc = input[xcc + (nx / 2) * (ny / 2 - 1 + (ny / 2) * zcc)],
					 urbc = input[xnc + (nx / 2) * (ny / 2 - 1 + (ny / 2) * zcc)],
					 ulfc = input[xcc + (nx / 2) * (ny / 2 - 1 + (ny / 2) * znc)],
					 urfc = input[xnc + (nx / 2) * (ny / 2 - 1 + (ny / 2) * znc)];
				auto &blbf = output[xcf + nx * ny * zcf],
				     &brbf = output[xnf + nx * ny * zcf],
				     &blff = output[xcf + nx * ny * znf],
				     &brff = output[xnf + nx * ny * znf],
				     &ulbf = output[xcf + ny * (ny - 1 + ny * zcf)],
				     &urbf = output[xnf + ny * (ny - 1 + ny * zcf)],
				     &ulff = output[xcf + ny * (ny - 1 + ny * znf)],
				     &urff = output[xnf + ny * (ny - 1 + ny * znf)];

				blbf = 9. / 16. * blbc + 3. / 16. * (brbc + blfc) + 1. / 16. * brfc;
				brbf = 9. / 16. * brbc + 3. / 16. * (blbc + brfc) + 1. / 16. * blfc;
				blff = 9. / 16. * blfc + 3. / 16. * (brfc + blbc) + 1. / 16. * brbc;
				brff = 9. / 16. * brfc + 3. / 16. * (blfc + brbc) + 1. / 16. * blbc;
				ulbf = 9. / 16. * ulbc + 3. / 16. * (urbc + ulfc) + 1. / 16. * urfc;
				urbf = 9. / 16. * urbc + 3. / 16. * (ulbc + urfc) + 1. / 16. * ulfc;
				ulff = 9. / 16. * ulfc + 3. / 16. * (urfc + ulbc) + 1. / 16. * urbc;
				urff = 9. / 16. * urfc + 3. / 16. * (ulfc + urbc) + 1. / 16. * ulbc;
			}
		}
	}

	template<class Geometry, unsigned Level>
	void
	solver<Geometry, Level>::restriction(
			dt::vector<nx * ny * nz / 8>& output,
			dt::vector<nx * ny * nz>& input)
	{
		for (auto xcc = 0u; xcc < nx / 2; ++xcc) {
			for (auto ycc = 0u; ycc < ny / 2; ++ycc) {
				for (auto zcc = 0u; zcc < nz / 2; ++zcc) {
					auto xcf = (2 * xcc + 1),
					     ycf = (2 * ycc + 1),
					     zcf = (2 * zcc + 1);
					auto& out = output[xcc + (nx / 2) * (ycc + (ny / 2) * zcc)];
					for (int i = 0; i < 64; ++i) {
						auto dx = (i % 4) - 2,
						     dy = (i % 16) / 4 - 2,
						     dz = i / 16 - 2;
						auto x = (xcf + dx + nx) % nx,
						     y = std::min(ycf + std::max(dy, -((int) ycf)), ny-1),
						     z = (zcf + dz + nz) % nz;
						auto weight = (4 - abs(2 * dx + 1)) *
						              (4 - abs(2 * dy + 1)) *
						              (4 - abs(2 * dz + 1));
						out += weight / 512. * input[x + nx * (y + ny * z)];
					}
				}
			}
		}
	}

	template<class Geometry, unsigned Level>
	void
	solver<Geometry, Level>::nested_iteration(dt::vector<nx * ny * nz>& input)
	{
		dt::vector<nx * ny * nz> copy(input), residual;
		dt::vector<nx * ny * nz / 8> coarse;

		restriction(coarse, input);
		coarser.nested_iteration(coarse);
		coarser.recurse(coarse);
		interpolation(input, coarse);
		
		for (auto j = 0u; j < nu2; ++j)
			alg::smooth(input, L, copy);
	}

	template<class Geometry, unsigned Level>
	void
	solver<Geometry, Level>::recurse(dt::vector<nx * ny * nz>& rhs)
	{
		dt::vector<nx * ny * nz> x, residualf;
		dt::vector<nx * ny * nz / 8> residualc;

		for (auto j = 0u; j < nu1; ++j)
			alg::smooth(x, L, rhs);
		residualf = rhs - L * x;
		restriction(residualc, residualf);
		coarser.recurse(residualc);
		interpolation(residualf, residualc);
		x += residualf;
		for (auto j = 0u; j < nu2; ++j)
			alg::smooth(x, L, rhs);
		rhs = x;
	}

	template<class Geometry, unsigned Level>
	unsigned int
	solver<Geometry, Level>::solve(dt::vector<nx * ny * nz>& rhs)
	{
		dt::vector<nx * ny * nz> x(rhs), residual, ones(1);
		nested_iteration(x);
		rhs -= (rhs * ones) / (nx * ny * nz) * ones;

		residual = rhs - L * x;
		auto count = 0u;
		while (abs(residual) > tol) {
			recurse(residual);
			x += residual;
			residual = rhs - L * x;
			++count;
		}
		rhs = x;

		return count;
	}

	template<class Geometry>
	class solver<Geometry, 0u> {
		public:
			static constexpr unsigned nx = Geometry::nx;
			static constexpr unsigned ny = Geometry::ny;
			static constexpr unsigned nz = Geometry::nz;
			static constexpr unsigned n = Geometry::n;

			typedef Geometry geometry_type;
			typedef laplacian<geometry_type> laplacian_type;
		private:
			double tol;
		protected:
			static laplacian_type L;
			void recurse(dt::vector<nx * ny * nz>& v) { solve(v); }
			void nested_iteration(dt::vector<nx * ny * nz>&) {};
		public:
			unsigned int solve(dt::vector<nx * ny * nz>&);
			solver(double tol, unsigned int, unsigned int) : tol(tol) {};

			friend class solver<typename geometry_type::fine>;
	};

	template<class Geometry> laplacian<Geometry> solver<Geometry, 0u>::L;

	template<class Geometry>
	unsigned int
	solver<Geometry, 0u>::solve(dt::vector<nx * ny * nz>& input)
	{
		dt::vector<nx * ny * nz> ones(1);
		input -= (ones * input) / (nx * ny * nz) * ones;
		return alg::conjugate_gradient(input, L, input, tol);
	}

}

#endif
