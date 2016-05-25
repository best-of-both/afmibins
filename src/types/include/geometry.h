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

#ifndef TYPES_GEOMETRY_H
#define TYPES_GEOMETRY_H

#include "typedefs.h"
#include "point.h"

namespace types {

	class geometry {
		private:
			static unsigned int compute_levels(size_type);
		public:
			static constexpr unsigned int factor = 2;
			const size_type width, height, depth, n;
			const size_type nx, ny, nz;

			point fit_to_box(point) const;
			bool fits_to_box(point) const;

			geometry(size_type w, size_type h, size_type d, size_type n) :
				width(w), height(h), depth(d), n(n), nx(w * n), ny(h * n), nz(d * n) {}
	};

}

#endif
