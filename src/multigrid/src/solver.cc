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
#include "types/geometry.h"
#include "types/vector.h"
#include "algo/smooth.h"
#include "algo/cg.h"
#include "solver.h"

namespace multigrid {

	using types::size_type;
	using types::index_type;
	using types::vector;

	solver::solver(double tol, unsigned int nu1, unsigned int nu2,
			const types::geometry& geometry) :
		m_top_solver(construct_solver(tol, nu1, nu2, geometry)) {}

	base_solver&
	solver::construct_solver(double tol, unsigned int nu1, unsigned int nu2,
			const types::geometry& g)
	{
		size_type n = g.n;
		unsigned int depth = 0;
		while (n % g.factor == 0) {
			n /= g.factor;
			++depth;
		}
		const types::geometry ng(g.width, g.height, g.depth, n);
		base_solver* current = new exact_solver(tol, nu1, nu2, ng);
		while (depth --> 0) {
			n *= g.factor;
			const types::geometry ng(g.width, g.height, g.depth, n);
			current = new iterative_solver(*current, ng);
		}
		return *current;
	}

	base_solver::base_solver(double tol, unsigned int nu1, unsigned int nu2,
			const types::geometry& geometry) :
		m_tol(tol), m_nu1(nu1), m_nu2(nu2), m_geometry(geometry), m_laplacian(m_geometry) {}

	exact_solver::exact_solver(double tol, unsigned int nu1, unsigned int nu2,
			const types::geometry& geometry) :
		base_solver(tol, nu1, nu2, geometry) {}

	unsigned int
	exact_solver::solve(vector& vec) const
	{
		const size_type size = m_geometry.nx * m_geometry.ny * m_geometry.nz;
		vector ones(size, 1);
		vec -= ((double) (vec * ones) / size) * ones;
		return algo::conjugate_gradient(vec, m_laplacian, vec, m_tol);
	}

	iterative_solver::iterative_solver(base_solver& child,
			const types::geometry& geometry) :
		base_solver(child.m_tol, child.m_nu1, child.m_nu2, geometry), m_child(child) {}

	void
	iterative_solver::interpolation(vector& out, vector& in) const
	{
		const unsigned int factor = m_geometry.factor;
		const size_type fnx = m_geometry.nx, fny = m_geometry.ny, fnz = m_geometry.nz;
		const size_type cnx = fnx / factor, cny = fny / factor, cnz = fnz / factor;

		for (index_type cxc = 0; cxc < cnx; ++cxc) {
			for (index_type cyc = 0; cyc < cny - 1; ++cyc) {
				for (index_type czc = 0; czc < cnz; ++czc) {
					const index_type cxn = (cxc + 1) % cnx,
					                 cyn = cyc + 1,
					                 czn = (czc + 1) % cnz;
					const index_type fxc = factor * cxc + 1,
					                 fyc = factor * cyc + 1,
					                 fzc = factor * czc + 1;
					const index_type fxn = (fxc + 1) % fnx,
					                 fyn = fyc + 1,
					                 fzn = (fzc + 1) % fnz;
					const double cmmm = in[cxc + cnx * (cyc + cny * czc)],
					             cpmm = in[cxn + cnx * (cyc + cny * czc)],
					             cmpm = in[cxc + cnx * (cyn + cny * czc)],
					             cppm = in[cxn + cnx * (cyn + cny * czc)],
					             cmmp = in[cxc + cnx * (cyc + cny * czn)],
					             cpmp = in[cxn + cnx * (cyc + cny * czn)],
					             cmpp = in[cxc + cnx * (cyn + cny * czn)],
					             cppp = in[cxn + cnx * (cyn + cny * czn)];
					double &fmmm = out[fxc + fnx * (fyc + fny * fzc)],
					       &fpmm = out[fxn + fnx * (fyc + fny * fzc)],
					       &fmpm = out[fxc + fnx * (fyn + fny * fzc)],
					       &fppm = out[fxn + fnx * (fyn + fny * fzc)],
					       &fmmp = out[fxc + fnx * (fyc + fny * fzn)],
					       &fpmp = out[fxn + fnx * (fyc + fny * fzn)],
					       &fmpp = out[fxc + fnx * (fyn + fny * fzn)],
					       &fppp = out[fxn + fnx * (fyn + fny * fzn)];

					fmmm = 27. / 64. * cmmm + 9. / 64. * (cpmm + cmpm + cmmp)
					     + 3. / 64. * (cppm + cpmp + cmpp) + 1. / 64. * cppp;
					fpmm = 27. / 64. * cpmm + 9. / 64. * (cmmm + cppm + cpmp)
					     + 3. / 64. * (cmpm + cmmp + cppp) + 1. / 64. * cmpp;
					fmpm = 27. / 64. * cmpm + 9. / 64. * (cppm + cmmm + cmpp)
					     + 3. / 64. * (cpmm + cmmp + cppp) + 1. / 64. * cpmp;
					fppm = 27. / 64. * cppm + 9. / 64. * (cpmm + cmpm + cppp)
					     + 3. / 64. * (cmmm + cpmp + cmpp) + 1. / 64. * cmmp;
					fmmp = 27. / 64. * cmmp + 9. / 64. * (cmmm + cmpp + cpmp)
					     + 3. / 64. * (cmpm + cpmp + cppp) + 1. / 64. * cppm;
					fpmp = 27. / 64. * cpmp + 9. / 64. * (cmmp + cppp + cpmm)
					     + 3. / 64. * (cmmm + cmpp + cppm) + 1. / 64. * cmpm;
					fmpp = 27. / 64. * cmpp + 9. / 64. * (cmmp + cmpm + cppp)
					     + 3. / 64. * (cmmm + cpmp + cppm) + 1. / 64. * cpmm;
					fppp = 27. / 64. * cppp + 9. / 64. * (cmpp + cpmp + cppm)
					     + 3. / 64. * (cmmp + cmpm + cpmm) + 1. / 64. * cmmm;
				}
			}

			for (index_type czc = 0; czc < cnz; ++czc) {
				const index_type cxn = (cxc + 1) % cnx,
				                 czn = (czc + 1) % cnz;
				const index_type fxc = factor * cxc + 1,
				                 fzc = factor * czc + 1;
				const index_type fxn = (fxc + 1) % fnx,
				                 fzn = (fzc + 1) % fnz;
				const double cmmm = in[cxc + cnx * cny * czc],
				             cpmm = in[cxn + cnx * cny * czc],
				             cmmp = in[cxc + cnx * cny * czn],
				             cpmp = in[cxn + cnx * cny * czn],
				             cmpm = in[cxc + cnx * (cny - 1 + cny * czc)],
				             cppm = in[cxn + cnx * (cny - 1 + cny * czc)],
				             cmpp = in[cxc + cnx * (cny - 1 + cny * czn)],
				             cppp = in[cxn + cnx * (cny - 1 + cny * czn)];
				double &fmmm = out[fxc + fnx * fny * fzc],
				       &fpmm = out[fxn + fnx * fny * fzc],
				       &fmmp = out[fxc + fnx * fny * fzn],
				       &fpmp = out[fxn + fnx * fny * fzn],
				       &fmpm = out[fxc + fnx * (fny - 1 + fny * fzc)],
				       &fppm = out[fxn + fnx * (fny - 1 + fny * fzc)],
				       &fmpp = out[fxc + fnx * (fny - 1 + fny * fzn)],
				       &fppp = out[fxn + fnx * (fny - 1 + fny * fzn)];

				fmmm = 9. / 16. * cmmm + 3. / 16. * (cmmp + cpmm) + 1. / 16. * cpmp;
				fpmm = 9. / 16. * cpmm + 3. / 16. * (cmmm + cpmp) + 1. / 16. * cmmp;
				fmmp = 9. / 16. * cmmp + 3. / 16. * (cmmm + cpmp) + 1. / 16. * cpmm;
				fpmp = 9. / 16. * cpmp + 3. / 16. * (cmmp + cpmm) + 1. / 16. * cmmm;
				fmpm = 9. / 16. * cmpm + 3. / 16. * (cmpp + cppm) + 1. / 16. * cppp;
				fppm = 9. / 16. * cppm + 3. / 16. * (cmpm + cppp) + 1. / 16. * cmpp;
				fmpp = 9. / 16. * cmpp + 3. / 16. * (cmpm + cppp) + 1. / 16. * cppm;
				fppp = 9. / 16. * cppp + 3. / 16. * (cmpp + cppm) + 1. / 16. * cmpm;
			}
		}
	}

	void
	iterative_solver::restriction(vector& out, vector& in) const
	{
		const unsigned int factor = m_geometry.factor;
		const size_type fnx = m_geometry.nx, fny = m_geometry.ny, fnz = m_geometry.nz;
		const size_type cnx = fnx / factor, cny = fny / factor, cnz = fnz / factor;

		for (index_type cxc = 0; cxc < cnx; ++cxc) {
			for (index_type cyc = 0; cyc < cny; ++cyc) {
				for (index_type czc = 0; czc < cnz; ++czc) {
					index_type fxc = factor * cxc + 1,
					           fyc = factor * cyc + 1,
					           fzc = factor * czc + 1;
					double& v = out[cxc + cnx * (cyc + cny * czc)];
					for (int i = 0; i < 64; ++i) {
						int dx = (i % 4) - 2,
						    dy = (i % 16) / 4 - 2,
							dz = i / 16 - 2;
						index_type x = (fxc + dx + fnx) % fnx,
						           y = std::min(std::max(fyc + dy, (size_type) 0), fny - 1),
						           z = (fzc + dz + fnz) % fnz;
						const double weight = (4 - abs(2 * dx + 1)) *
						                      (4 - abs(2 * dy + 1)) *
						                      (4 - abs(2 * dz + 1));
						v += weight / 512. * in[x + fnx * (y + fny * z)];
					}
				}
			}
		}
	}

	void
	iterative_solver::nested_iteration(vector& in) const
	{
		const unsigned int factor = m_geometry.factor;
		const size_type size = m_geometry.nx * m_geometry.ny * m_geometry.nz;
		vector copy(size), residual(size);
		vector coarse(size / (factor * factor * factor));

		copy = in;
		restriction(coarse, in);
		m_child.nested_iteration(coarse);
		m_child.recurse(coarse);
		interpolation(in, coarse);
		for (unsigned int i = 0; i < m_nu2; ++i)
			algo::smooth(in, m_laplacian, copy);
	}

	void
	iterative_solver::recurse(vector& rhs) const
	{
		const unsigned int factor = m_geometry.factor;
		const size_type size = m_geometry.nx * m_geometry.ny * m_geometry.nz;
		vector x(size), rfine(size);
		vector rcoarse(size / (factor * factor * factor));

		for (unsigned int i = 0; i < m_nu1; ++i)
			algo::smooth(x, m_laplacian, rhs);
		rfine = rhs - m_laplacian * x;
		restriction(rcoarse, rfine);
		m_child.recurse(rcoarse);
		interpolation(rfine, rcoarse);
		x += rfine;
		for (unsigned int i = 0; i < m_nu2; ++i)
			algo::smooth(x, m_laplacian, rhs);
		rhs = x;
	}

	unsigned int
	iterative_solver::solve(vector& rhs) const
	{
		const size_type size = m_geometry.nx * m_geometry.ny * m_geometry.nz;
		vector x(size), residual(size), ones(size, 1);

		x = rhs;
		rhs -= (((double) (ones * rhs)) / size) * ones;
		nested_iteration(x);

		residual = rhs - m_laplacian * x;
		unsigned int count = 0;
		while (abs(residual) > m_tol) {
			recurse(residual);
			x += residual;
			residual = rhs - m_laplacian * x;
			++count;
		}
		rhs = x;

		return count;
	}

}
