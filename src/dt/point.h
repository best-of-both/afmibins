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

#ifndef DT_POINT_H
#define DT_POINT_H

#include <ostream>

namespace dt {

	class point {
		public:
			typedef double value_type;
			typedef value_type& reference;
			typedef const value_type& const_reference;
			typedef std::size_t size_type;

		private:
			value_type values[3];
			size_type index;
		public:
			point(value_type x, value_type y, value_type z, const size_type index)
				: values{x, y, z}, index(index) {}
			point() : values{0.0, 0.0, 0.0}, index(0) {}

			reference operator[](size_type index) { return values[index]; }
			const_reference operator[](size_type index) const { return values[index]; }
			size_type get_index() const { return index; }

	};

	std::ostream&
	operator<<(std::ostream& out, point p)
	{
		return out << p.get_index() << " (" << p[0] << ", " << p[1] << ", " << p[2] << ")";
	}

}

#endif
