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

#ifndef DT_GEOMETRY_H
#define DT_GEOMETRY_H

#include "vector.h"

namespace dt {

	template<unsigned Level, unsigned Width, unsigned Height, unsigned Depth>
	class geometry {
		public:
			typedef geometry<Level-1, Width, Height, Depth> coarse;
			typedef geometry<Level+1, Width, Height, Depth> fine;
			static constexpr unsigned factor = 2;
			static constexpr unsigned n = coarse::n * factor;
			static constexpr unsigned nx = Width * n;
			static constexpr unsigned ny = Height * n;
			static constexpr unsigned nz = Depth * n;
			static constexpr unsigned size = nx * ny * nz;
			static constexpr unsigned width = Width;
			static constexpr unsigned height = Height;
			static constexpr unsigned depth = Depth;
			static constexpr unsigned level = Level;
	};

	template<unsigned Width, unsigned Height, unsigned Depth>
	class geometry<0u, Width, Height, Depth> {
		public:
			typedef geometry<1u, Width, Height, Depth> fine;
			static constexpr unsigned factor = 2;
			static constexpr unsigned n = 1;
			static constexpr unsigned nx = Width;
			static constexpr unsigned ny = Height;
			static constexpr unsigned nz = Depth;
			static constexpr unsigned size = nx * ny * nz;
			static constexpr unsigned width = Width;
			static constexpr unsigned height = Height;
			static constexpr unsigned depth = Depth;
			static constexpr unsigned level = 0u;
	};

}

#endif
