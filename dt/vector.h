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

#ifndef DT_VECTOR_H
#define DT_VECTOR_H

#include <vector>
#include <cstring>
#include <cmath>
#include <ostream>
#include <cassert>

namespace dt {

	template<unsigned int Size>
	class vector {
		private:
			typedef std::vector<double> vector_type;
		public:
			typedef typename vector_type::value_type value_type;
			typedef typename vector_type::iterator iterator;
			typedef typename vector_type::const_iterator const_iterator;
			typedef typename vector_type::reference reference;
			typedef typename vector_type::const_reference const_reference;
			typedef typename vector_type::pointer pointer;
			typedef typename vector_type::const_pointer const_pointer;
			typedef typename vector_type::size_type size_type;
		protected:
			vector_type values;
		public:
			value_type operator*(vector<Size>) const;
			vector<Size>& operator+=(vector<Size>);
			vector<Size>& operator-=(vector<Size>);
			vector<Size>& operator*=(double);
			vector<Size>& operator/=(double);
			vector<Size>& operator=(std::initializer_list<value_type> v)
			{
				assert(v.size() == Size);
				values = v;
			}
			virtual reference operator[](unsigned int i) { return values[i]; }
			virtual const_reference operator[](unsigned int i) const { return values[i]; }
			size_type size() const { return Size; };
			iterator begin() { return values.begin(); }
			iterator end() { return values.end(); }
			const_iterator cbegin() const { return values.cbegin(); }
			const_iterator cend() const { return values.cend(); }

			vector(const vector<Size>& vec)
				: values(vec.values) {};
			vector(double value = 0.0)
				: values(Size, value) {};
			vector(std::initializer_list<value_type> values)
				: values(values) { assert(values.size() == Size); };
	};

	template<unsigned int Size>
	double
	abs(vector<Size>& v)
	{
		double mod = 0.0;
		for (const auto& e: v)
			mod += e * e;
		return sqrt(mod) / sqrt(Size);
	}

	template<unsigned int Size>
	auto
	vector<Size>::operator*(vector<Size> other) const ->
		value_type
	{
		value_type result = 0.0;
		for (auto i = 0u; i < Size; ++i)
			result += values[i] * other.values[i];
		return result;
	}

	template<unsigned int Size>
	vector<Size>&
	vector<Size>::operator*=(double c)
	{
		for (auto& v: *this)
			v *= c;
		return *this;
	}

	template<unsigned int Size>
	vector<Size>&
	vector<Size>::operator/=(double c) {
		for (auto& v: *this)
			v /= c;
		return *this;
	}

	template<unsigned int Size>
	vector<Size>&
	vector<Size>::operator+=(vector<Size> other) {
		for (auto i = 0u; i < Size; ++i)
			values[i] += other[i];
		return *this;
	}

	template<unsigned int Size>
	vector<Size>&
	vector<Size>::operator-=(vector<Size> other) {
		for (auto i = 0u; i < Size; ++i)
			values[i] -= other[i];
		return *this;
	}

	template<unsigned int Size>
	vector<Size>
	operator+(vector<Size> left, vector<Size> right) {
		return left += right;
	}

	template<unsigned int Size>
	vector<Size>
	operator-(vector<Size> left, vector<Size> right) {
		return left -= right;
	}

	template<unsigned int Size>
	vector<Size>
	operator*(double c, vector<Size> vec) {
		return vec *= c;
	}

	template<unsigned int Size>
	vector<Size>
	operator/(vector<Size> vec, double c) {
		return vec /= c;
	}

	template<unsigned int Size>
	std::ostream&
	operator<<(std::ostream& out, vector<Size> vec) {
		bool p = false;
		out << "[";
		for (auto& v: vec) {
			if (p) out << ", ";
			p = true;
			out << v;
		}
		out << "]'";
		return out;
	}
}

#endif
