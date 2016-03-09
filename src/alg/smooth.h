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

#ifndef ALG_SMOOTH_H
#define ALG_SMOOTH_H

namespace dt {
	template<unsigned int>
	class matrix;

	template<unsigned int>
	class vector;
}

namespace alg {

	template<unsigned int Size>
	void
	smooth(dt::vector<Size>& x, const dt::matrix<Size>& A, const dt::vector<Size>& b)
	{
		for (std::size_t row = 0; row < Size; ++row) {
			typename dt::vector<Size>::value_type value = 0.0;
			typename dt::matrix<Size>::value_type diagonal = 0.0;
			for (auto i = A.row_start[row]; i < A.row_start[row+1]; ++i) {
				if (A.column[i] == row) {
					diagonal = A.values[i];
					continue;
				}
				value += A.values[i] * x[A.column[i]];
			}
			x[row] = (b[row] - value) / diagonal;
		}
	}

}

#endif
