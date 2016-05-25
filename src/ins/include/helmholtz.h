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

#ifndef INS_HELMHOLTZ_H
#define INS_HELMHOLTZ_H

#include "types/typedefs.h"
#include "types/geometry.h"
#include "types/vector.h"
#include "types/matrix.h"

namespace ins {

	using types::size_type;
	using types::vector;

	class helmholtz_periodic : public types::sparse_square_matrix {
		protected:
			const types::geometry& m_geometry;
			const double m_scale;
		public:
			vector boundary(vector&, vector&) const;

			helmholtz_periodic(const types::geometry&, double);
	};

	class helmholtz_dirichlet : public types::sparse_square_matrix {
		protected:
			const types::geometry& m_geometry;
			const double m_scale;
		public:
			vector boundary(vector&, vector&) const;

			helmholtz_dirichlet(const types::geometry&, double);
	};

}

#endif
