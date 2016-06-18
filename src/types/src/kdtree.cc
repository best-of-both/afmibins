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

#include "point.h"
#include "kdtree.h"

namespace types {

	std::vector<point>
	kdtree_base::sort(iterator start, iterator end, bool (*cmp)(point, point))
	{
		std::sort(start, end, cmp);
		return std::vector<point>(start, end);
	}

	kdtree_x::kdtree_x(std::vector<point> l, const geometry& g) :
		kdtree_base(g.fit_to_box(l[l.size()/2]),
			l.size() > 1 ? new tree_t(sort(l.begin(), l.begin() + l.size()/2, tree_t::comp), g) : nullptr,
			l.size() > 2 ? new tree_t(sort(l.begin() + l.size()/2 + 1, l.end(), tree_t::comp), g) : nullptr,
			g.width, g.n) {}

	std::vector<point>
	kdtree_x::search_nearby(point coords) const
	{
		return search_nearby_periodic(coords, *this);
	}

	kdtree_y::kdtree_y(std::vector<point> l, const geometry& g) :
		kdtree_base(l[l.size() / 2],
			l.size() > 1 ? new tree_t(sort(l.begin(), l.begin() + l.size()/2, tree_t::comp), g) : nullptr,
			l.size() > 2 ? new tree_t(sort(l.begin() + l.size()/2 + 1, l.end(), tree_t::comp), g) : nullptr,
			g.height, g.n) {}

	std::vector<point>
	kdtree_y::search_nearby(point coords) const
	{
		return search_nearby_nonperiodic(coords, *this);
	}

	kdtree_z::kdtree_z(std::vector<point> l, const geometry& g) :
		kdtree_base(l[l.size() / 2],
			l.size() > 1 ? new tree_t(sort(l.begin(), l.begin() + l.size()/2, tree_t::comp), g) : nullptr,
			l.size() > 2 ? new tree_t(sort(l.begin() + l.size()/2 + 1, l.end(), tree_t::comp), g) : nullptr,
			g.depth, g.n) {}

	std::vector<point>
	kdtree_z::search_nearby(point coords) const
	{
		return search_nearby_periodic(coords, *this);
	}

	kdtree::kdtree(std::vector<point> l, const geometry& g) :
		m_tree(new tree_t(kdtree_base::sort(l.begin(), l.end(), tree_t::comp), g)) {}

	std::vector<point>
	kdtree::search_nearby(point coords) const
	{
		return m_tree->search_nearby(coords);
	}

	std::vector<point>
	search_nearby_nonperiodic(point coords, const kdtree_base& tree)
	{
		std::vector<point> results;
		const point& ref = tree.m_location;
		const kdtree_base *left = tree.m_left, *right = tree.m_right;
		const size_type &n = tree.m_n;

		const double re = tree.entry(ref), ce = tree.entry(coords);
		const double de = ce - re;
		const double threshold = 2. / n;

		if (de <= -threshold) {
			if (left != nullptr)
				for (point& result: left->search_nearby(coords))
					results.push_back(result);
		} else if (de >= threshold) {
			if (right != nullptr)
				for (point& result: right->search_nearby(coords))
					results.push_back(result);
		} else {
			if (right != nullptr)
				for (point& result: right->search_nearby(coords))
					results.push_back(result);
			if (left != nullptr)
				for (point& result: left->search_nearby(coords))
					results.push_back(result);

			point diff = ref - coords;
			if (std::abs(diff.x) >= threshold ||
					std::abs(diff.y) >= threshold ||
			 		std::abs(diff.z) >= threshold) {
				return results;
			}
			results.push_back(ref);
		}
		return results;
	}

	std::vector<point>
	search_nearby_periodic(point coords, const kdtree_base& tree)
	{
		std::vector<point> results = search_nearby_nonperiodic(coords, tree);
		const point& ref = tree.m_location;
		const kdtree_base *left = tree.m_left, *right = tree.m_right;
		const size_type &length = tree.m_length;
		const size_type &n = tree.m_n;

		const double re = tree.entry(ref);
		double& ce = tree.entry(coords);
		const double de = ce - re;
		const double threshold = 2. / n;

		if (de < -threshold && ce < threshold && right != nullptr) {
			ce += length;
			for (point result: right->search_nearby(coords)) {
				tree.entry(result) -= length;
				results.push_back(result);
			}
		}
		else if (de > threshold && ce > length - threshold && left != nullptr) {
			ce -= length;
			for (point result: left->search_nearby(coords)) {
				tree.entry(result) += length;
				results.push_back(result);
			}
		}

		return results;
	}

}
