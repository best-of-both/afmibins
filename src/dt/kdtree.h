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

#ifndef DT_KDTREE_H
#define DT_KDTREE_H

#include <algorithm>
#include <iterator>
#include <vector>
#include <utility>
#include <cmath>
#include <memory>

#include "point.h"
#include "geometry.h"

namespace dt {

	template<class, unsigned int = 0> class kdtree;

	template<class Tree>
	std::vector<typename Tree::point_type>
	search_nearby_periodic(typename Tree::point_type coords, const Tree& tree)
	{
		std::vector<typename Tree::point_type> results;
		auto& ref = tree.location;
		auto *left = tree.left, *right = tree.right;
		auto length = tree.length;
		auto axis = tree.axis;
		auto n = tree.n;


		if (coords[axis] + 2. / n < ref[axis]) {
			if (left != nullptr) {
				for (auto& result: left->search_nearby(coords))
					results.push_back(result);
			}

			if (coords[axis] < 2. / n) {
				if (right != nullptr) {
					coords[axis] += length;
					for (auto result: right->search_nearby(coords)) {
						result[axis] -= length;
						results.push_back(result);
					}
				}
			}
		}
		else if (coords[axis] - 2. / n > ref[axis]) {
			if (right != nullptr) {
				for (auto& result: right->search_nearby(coords))
					results.push_back(result);
			}

			if (coords[axis] > length - 2. / n) {
				if (left != nullptr) {
					coords[axis] -= length;
					for (auto& result: left->search_nearby(coords)) {
						result[axis] += length;
						results.push_back(result);
					}
				}
			}
		} else {
			if (right != nullptr) {
				for (auto& result: right->search_nearby(coords))
					results.push_back(result);
			}
			if (left != nullptr) {
				for (auto& result: left->search_nearby(coords))
					results.push_back(result);
			}

			for (std::size_t i = 0; i < 3; ++i) {
				auto local = ref[i];
				auto remote = coords[i];
				if (std::abs(local - remote) > 2. / n)
					return results;
			}

			results.push_back(ref);
		}
		return results;
	}

	template<class Tree>
	std::vector<typename Tree::point_type>
	search_nearby_nonperiodic(typename Tree::point_type coords, const Tree& tree)
	{
		std::vector<typename Tree::point_type> results;
		auto& ref = tree.location;
		auto *left = tree.left, *right = tree.right;
		auto axis = tree.axis;
		auto n = tree.n;


		if (coords[axis] + 2. / n < ref[axis]) {
			if (left != nullptr) {
				for (auto& result: left->search_nearby(coords))
					results.push_back(result);
			}
		}
		else if (coords[axis] - 2. / n > ref[axis]) {
			if (right != nullptr) {
				for (auto& result: right->search_nearby(coords))
					results.push_back(result);
			}
		} else {
			if (right != nullptr) {
				for (auto& result: right->search_nearby(coords))
					results.push_back(result);
			}
			if (left != nullptr) {
				for (auto& result: left->search_nearby(coords))
					results.push_back(result);
			}

			for (std::size_t i = 0; i < 3; ++i) {
				auto local = ref[i];
				auto remote = coords[i];
				if (std::abs(local - remote) > 2. / n)
					return results;
			}

			results.push_back(ref);
		}
		return results;
	}

	template<class Geometry>
	class kdtree<Geometry, 0> {
		public:
			typedef point point_type;
			typedef kdtree<Geometry, 1>* tree_pointer;
		private:
			class __sorter {
				public:
					bool operator()(point_type left, point_type right) { return left[0] < right[0]; }
			} sorter;
		protected:
			point_type location;
			static constexpr unsigned int n = Geometry::n;
			static constexpr unsigned int length = Geometry::width;
			static constexpr unsigned int axis = 0;
			tree_pointer left, right;
		public:
			kdtree(std::vector<point_type> list)
				: left(nullptr), right(nullptr)
			{
				auto size = list.size();
				auto index = size / 2;

				std::sort(list.begin(), list.end(), sorter);
				location = list[index];

				if (index > 0) {
					std::vector<point_type> points(list.begin(), list.begin() + index);
					left = new kdtree<Geometry, 1>(points);
				}
				if (index + 1 < size) {
					std::vector<point_type> points(list.begin() + index + 1, list.end());
					right = new kdtree<Geometry, 1>(points);
				}
			}
			~kdtree() { delete left; delete right; }

			std::vector<point_type> search_nearby(point_type coords) const
			{
				return search_nearby_periodic(coords, *this);
			}

			friend std::vector<point_type> search_nearby_periodic<kdtree<Geometry, 0>>(point_type, const kdtree<Geometry, 0> &);
	};

	template<class Geometry>
	class kdtree<Geometry, 1> {
		public:
			typedef point point_type;
			typedef kdtree<Geometry, 2>* tree_pointer;
		private:
			class __sorter {
				public:
					bool operator()(point_type left, point_type right) { return left[1] < right[1]; }
			} sorter;
		protected:
			point_type location;
			static constexpr unsigned int n = Geometry::n;
			static constexpr unsigned int length = Geometry::height;
			static constexpr unsigned int axis = 1;
			kdtree<Geometry, 2> *left, *right;

		public:
			kdtree(std::vector<point_type> list)
				: left(nullptr), right(nullptr)
			{
				auto size = list.size();
				auto index = size / 2;

				std::sort(list.begin(), list.end(), sorter);
				location = list[index];

				if (index > 0) {
					std::vector<point_type> points(list.begin(), list.begin() + index);
					left = new kdtree<Geometry, 2>(points);
				}
				if (index + 1 < size) {
					std::vector<point_type> points(list.begin() + index + 1, list.end());
					right = new kdtree<Geometry, 2>(points);
				}
			}
			~kdtree() { delete left; delete right; }

			std::vector<point_type> search_nearby(point_type coords) const
			{
				return search_nearby_nonperiodic(coords, *this);
			}

			friend std::vector<point_type> search_nearby_nonperiodic<kdtree<Geometry, 1>>(point_type, const kdtree<Geometry, 1>&);
	};

	template<class Geometry>
	class kdtree<Geometry, 2> {
		public:
			typedef point point_type;
			typedef kdtree<Geometry, 1>* tree_pointer;
		private:
			class __sorter {
				public:
					bool operator()(point_type left, point_type right) { return left[2] < right[2]; }
			} sorter;
		protected:
			point_type location;
			static constexpr unsigned int n = Geometry::n;
			static constexpr unsigned int length = Geometry::depth;
			static constexpr unsigned int axis = 2;
			kdtree<Geometry, 0> *left, *right;

		public:
			kdtree(std::vector<point_type> list)
				: left(nullptr), right(nullptr)
			{
				auto size = list.size();
				auto index = size / 2;

				std::sort(list.begin(), list.end(), sorter);
				location = list[index];

				if (index > 0) {
					std::vector<point_type> points(list.begin(), list.begin() + index);
					left = new kdtree<Geometry, 0>(points);
				}
				if (index + 1 < size) {
					std::vector<point_type> points(list.begin() + index + 1, list.end());
					right = new kdtree<Geometry, 0>(points);
				}
			}
			~kdtree() { delete left; delete right; }

			std::vector<point_type> search_nearby(point_type coords) const
			{
				return search_nearby_periodic(coords, *this);
			}

			friend std::vector<point_type> search_nearby_periodic<kdtree<Geometry, 2>>(point_type, const kdtree<Geometry, 2>&);
	};
}

#endif
