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

#ifndef TYPES_POINT_H
#define TYPES_POINT_H

#include <ostream>

#include "typedefs.h"

namespace types {

	class point {
		public:
			double x, y;
			index_type index;
			point& operator+=(point);
			point& operator-=(point);
			point& operator*=(double);
			point& operator/=(double);

			point(double x, double y, index_type index = 0) :
				x(x), y(y), index(index) {}
	};

	point operator+(point, point);
	point operator-(point, point);
	point operator*(double, point);
	point operator*(point, double);
	point operator/(point, double);

	double abs(point);
	std::ostream& operator<<(std::ostream&, point);

}

#endif
