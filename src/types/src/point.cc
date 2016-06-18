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

#include <cmath>
#include <ostream>

#include "point.h"

namespace types {

	point&
	point::operator+=(point other)
	{
		x += other.x;
		y += other.y;
		z += other.z;
		return *this;
	}

	point&
	point::operator-=(point other)
	{
		x -= other.x;
		y -= other.y;
		z -= other.z;
		return *this;
	}

	point&
	point::operator*=(double other)
	{
		x *= other;
		y *= other;
		z *= other;
		return *this;
	}

	point&
	point::operator/=(double other)
	{
		x /= other;
		y /= other;
		z /= other;
		return *this;
	}

	point
	operator+(point left, point right)
	{
		return left += right;
	}

	point
	operator-(point left, point right)
	{
		return left -= right;
	}

	point
	operator*(double left, point right)
	{
		return right *= left;
	}

	point 
	operator*(point left, double right)
	{
		return left *= right;
	}

	point
	operator/(point left, double right)
	{
		return left /= right;
	}

	double
	abs(point p)
	{
		return sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
	}

	std::ostream&
	operator<<(std::ostream& out, point p)
	{
		return out << "(" << p.x << ", " << p.y << ", " << p.z << ")";
	}
}

using types::operator+;
using types::operator-;
using types::operator*;
using types::operator/;
