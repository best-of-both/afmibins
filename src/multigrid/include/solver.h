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

#ifndef MULTIGRID_SOLVER_H
#define MULTIGRID_SOLVER_H

#include <iostream>

#include "types/typedefs.h"
#include "types/geometry.h"
#include "types/vector.h"
#include "laplacian.h"

#define FACTOR 2

namespace multigrid {

	using types::size_type;
	using types::vector;

	class solver;
	class exact_solver;
	class iterative_solver;

	class base_solver {
		protected:
			const double m_tol;
			const unsigned int m_nu1, m_nu2;
			const types::geometry m_geometry;
			laplacian m_laplacian;

			virtual void nested_iteration(vector&) const = 0;
			virtual void recurse(vector&) const = 0;
			virtual unsigned int solve(vector&) const = 0;

			base_solver(double, unsigned int, unsigned int,
					const types::geometry&);
			base_solver(base_solver&, size_type n);
			virtual ~base_solver() {}

		friend class solver;
		friend class exact_solver;
		friend class iterative_solver;
	};

	class exact_solver : public base_solver {
		protected:
			using base_solver::m_laplacian;

			void nested_iteration(vector&) const {};
			void recurse(vector& v) const { solve(v); }
			unsigned int solve(vector&) const;

			exact_solver(double, unsigned int, unsigned int,
					const types::geometry&);

		friend class solver;
		friend class iterative_solver;
	};

	class iterative_solver : public base_solver {
		private:
			const base_solver& m_child;

			void restriction(vector&, vector&) const;
			void interpolation(vector&, vector&) const;
		protected:
			using base_solver::m_laplacian;

			void nested_iteration(vector&) const;
			void recurse(vector&) const;
			unsigned int solve(vector&) const;

			iterative_solver(base_solver&,
					const types::geometry&);
			virtual ~iterative_solver() { delete &m_child; }

		friend class solver;
	};

	class solver {
		private:
			const base_solver& m_top_solver;
			static base_solver& construct_solver(double, unsigned int, unsigned int,
					const types::geometry&);
		public:
			void solve(vector& v) const { m_top_solver.solve(v); }

			solver(double, unsigned int, unsigned int, const types::geometry&);
			~solver() { delete &m_top_solver; }
	};
}

#endif
