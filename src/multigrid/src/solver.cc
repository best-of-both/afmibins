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
		const types::geometry ng(g.width, g.height, n);
		base_solver* current = new exact_solver(tol, nu1, nu2, ng);
		while (depth --> 0) {
			n *= g.factor;
			const types::geometry ng(g.width, g.height, n);
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
		const size_type size = m_geometry.nx * m_geometry.ny;
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
		const size_type fnx = m_geometry.nx, fny = m_geometry.ny;
		const size_type cnx = fnx / factor, cny = fny / factor;

		for (index_type cxc = 0; cxc < cnx; ++cxc) {
			const index_type cxn = (cxc + 1) % cnx;
			const index_type fxc = factor * cxc + 1;
			const index_type fxn = (fxc + 1) % fnx;
			for (index_type cyc = 0; cyc < cny - 1; ++cyc) {
				const index_type cyn = cyc + 1;
				const index_type fyc = factor * cyc + 1;
				const index_type fyn = fyc + 1;
				const double cmm = in[cxc + cnx * cyc],
							 cpm = in[cxn + cnx * cyc],
							 cmp = in[cxc + cnx * cyn],
							 cpp = in[cxn + cnx * cyn];
				double &fmm = out[fxc + fnx * fyc],
					   &fpm = out[fxn + fnx * fyc],
					   &fmp = out[fxc + fnx * fyn],
					   &fpp = out[fxn + fnx * fyn];

				fmm = 9. / 16. * cmm + 3. / 16. * (cpm + cmp) + 1. / 16. * cpp;
				fpm = 9. / 16. * cpm + 3. / 16. * (cmm + cpp) + 1. / 16. * cmp;
				fmp = 9. / 16. * cmp + 3. / 16. * (cpp + cmm) + 1. / 16. * cpm;
				fpp = 9. / 16. * cpp + 3. / 16. * (cpm + cmp) + 1. / 16. * cmm;
			}

			const double cmm = in[cxc],
						 cpm = in[cxn],
						 cmp = in[cxc + cnx * (cny - 1)],
						 cpp = in[cxn + cnx * (cny - 1)];
			double &fmm = out[fxc],
				   &fpm = out[fxn],
				   &fmp = out[fxc + fnx * (fny - 1)],
				   &fpp = out[fxn + fnx * (fny - 1)];

			fmm = 3. / 4. * cmm + 1. / 4. * cpm;
			fpm = 3. / 4. * cpm + 1. / 4. * cmm;
			fmp = 3. / 4. * cmp + 1. / 4. * cpp;
			fpp = 3. / 4. * cpp + 1. / 4. * cmp;
		}
	}

	void
	iterative_solver::restriction(vector& out, vector& in) const
	{
		const unsigned int factor = m_geometry.factor;
		const size_type fnx = m_geometry.nx, fny = m_geometry.ny;
		const size_type cnx = fnx / factor, cny = fny / factor;

		for (index_type cxc = 0; cxc < cnx; ++cxc) {
			index_type fxc = factor * cxc + 1;
			for (index_type cyc = 0; cyc < cny; ++cyc) {
				index_type fyc = factor * cyc + 1;
				double& v = out[cxc + cnx * cyc];
				for (int i = 0; i < 16; ++i) {
					int dx = i % 4 - 2,
					    dy = i / 4 - 2;
					index_type x = (fxc + dx + fnx) % fnx,
							   y = std::min(std::max(fyc + dy, (size_type) 0), fny - 1);
					const double weight = (4 - abs(2 * dx + 1)) *
					                      (4 - abs(2 * dy + 1));
					v += weight / 64. * in[x + fnx * y];
				}
			}
		}
	}

	void
	iterative_solver::nested_iteration(vector& in) const
	{
		const unsigned int factor = m_geometry.factor;
		const size_type size = m_geometry.nx * m_geometry.ny;
		vector copy(size), residual(size);
		vector coarse(size / (factor * factor));

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
		const size_type size = m_geometry.nx * m_geometry.ny;
		vector x(size), rfine(size);
		vector rcoarse(size / (factor * factor));

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
		const size_type size = m_geometry.nx * m_geometry.ny;
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
