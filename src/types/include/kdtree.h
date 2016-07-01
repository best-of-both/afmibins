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

#ifndef TYPES_KDTREE_H
#define TYPES_KDTREE_H

#include <algorithm>
#include <vector>
#include <cmath>

#include "geometry.h"
#include "point.h"

namespace types {

	class kdtree;
	class kdtree_base;
	class kdtree_x;
	class kdtree_y;

	std::vector<point> search_nearby_periodic(point, const kdtree_base&);
	std::vector<point> search_nearby_nonperiodic(point, const kdtree_base&);

	class kdtree_base {
		protected:
			typedef std::vector<point>::iterator iterator;

			const point m_location;
			const kdtree_base *m_left, *m_right;
			const size_type m_length, m_n;

			virtual double& entry(point&) const = 0;
			virtual const double& entry(const point&) const = 0;
			virtual std::vector<point> search_nearby(point) const = 0;
			static std::vector<point> sort(iterator, iterator, bool (*)(point, point));

			kdtree_base(point median, kdtree_base* left, kdtree_base* right,
					size_type length, size_type n) :
				m_location(median), m_left(left), m_right(right), m_length(length), m_n(n) {}
			virtual ~kdtree_base() { delete m_left; delete m_right; }
		friend std::vector<point> search_nearby_periodic(point, const kdtree_base&);
		friend std::vector<point> search_nearby_nonperiodic(point, const kdtree_base&);
		friend class kdtree;
	};

	class kdtree_x : public kdtree_base {
		private:
			typedef kdtree_y tree_t;
		protected:
			static bool comp(point l, point r) { return l.x < r.x; }
			double& entry(point& p) const { return p.x; }
			const double& entry(const point& p) const { return p.x; }
			std::vector<point> search_nearby(point) const;

			kdtree_x(std::vector<point>, const geometry&);

		friend class kdtree;
		friend class kdtree_y;
	};

	class kdtree_y : public kdtree_base {
		private:
			typedef kdtree_x tree_t;
		protected:
			static bool comp(point l, point r) { return l.y < r.y; }
			double& entry(point& p) const { return p.y; }
			const double& entry(const point& p) const { return p.y; }
			std::vector<point> search_nearby(point) const;

			kdtree_y(std::vector<point>, const geometry&);

		friend class kdtree_x;
	};

	class kdtree {
		private:
			typedef kdtree_x tree_t;
		protected:
			const kdtree_x* m_tree;
		public:
			std::vector<point> search_nearby(point) const;

			kdtree(std::vector<point>, const geometry&);
			~kdtree() { delete m_tree; }
	};
}

#endif
