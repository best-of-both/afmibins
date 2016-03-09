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

#ifndef HELMHOLTZ_H
#define HELMHOLTZ_H

#include "dt/matrix.h"

namespace ns {

    template<class Geometry, unsigned int Size>
    class helmholtz : public dt::matrix<Size> {
        protected:
            using dt::matrix<Size>::values;
            using dt::matrix<Size>::column;

            helmholtz() : dt::matrix<Size>() {
                values.reserve(Size * std::min(7u, Size));
                column.reserve(Size * std::min(7u, Size));
            };
    };

    template<class Geometry>
    class helmholtz_periodic : public helmholtz<Geometry, Geometry::nx * Geometry::ny * Geometry::nz> {
		public:
			static constexpr unsigned int nx = Geometry::nx;
			static constexpr unsigned int ny = Geometry::ny;
			static constexpr unsigned int nz = Geometry::nz;
			static constexpr unsigned int n = Geometry::n;
        private:
            double scale;
        protected:
            using dt::matrix<nx * ny * nz>::values;
            using dt::matrix<nx * ny * nz>::column;
            using dt::matrix<nx * ny * nz>::row_start;
        public:
            dt::vector<nx * ny * nz> boundary(
                    dt::vector<nx * nz>&,
                    dt::vector<nx * nz>&);
            helmholtz_periodic(double);
    };

    template<class Geometry>
    helmholtz_periodic<Geometry>::helmholtz_periodic(double k)
        : helmholtz<Geometry, nx * ny * nz> (), scale(k)
    {
        for (auto i = 0u; i < nx * ny * nz; ++ i) {
            row_start[i] = values.size();
            auto x = i % nx, y = (i % (nx * ny)) / nx, z = i / (nx * ny);
            const double offdiagonal = scale * n * n;
            double diagonal = 1. - 6. * offdiagonal,
                   left = offdiagonal,
                   right = offdiagonal,
                   top = offdiagonal,
                   bottom = offdiagonal,
                   front = offdiagonal,
                   back = offdiagonal;

            if (nx == 1) {
                diagonal += left + right;
                left = right = 0;
            } else if (nx == 2) {
                if (x == 0) {
                    right += left;
                    left = 0;
                } else {
                    left += right;
                    right = 0;
                }
            }

            if (ny == 1) {
                diagonal -= 3 * (top + bottom);
                top = bottom = 0;
            } else {
                if (y == 0) {
                    diagonal -= 2 * bottom;
                    top += offdiagonal / 3;
                    bottom = 0;
                }
                else if (y == ny - 1) {
                    diagonal -= 2 * top;
                    bottom += offdiagonal / 3;
                    top = 0;
                }
            }

            if (nz == 1) {
                diagonal += back + front;
                back = front = 0;
            } else if (nz == 2) {
                if (z == 0) {
                    front += back;
                    back = 0;
                } else {
                    back += front;
                    front = 0;
                }
            }

            if (front != 0 && z == nz - 1) {
                values.push_back(front);
                column.push_back(x + nx * y);
            }

            if (back != 0 && z > 0) {
                values.push_back(back);
                column.push_back(i - nx * ny);
            }

            if (bottom != 0 && y > 0) {
                values.push_back(bottom);
                column.push_back(i - nx);
            }

            if (right != 0 && x == nx - 1) {
                values.push_back(right);
                column.push_back(nx * (y + ny * z));
            }

            if (left != 0 && x > 0) {
                values.push_back(left);
                column.push_back(i - 1);
            }

            if (diagonal != 0) {
                values.push_back(diagonal);
                column.push_back(i);
            }

            if (right != 0 && x < nx - 1) {
                values.push_back(right);
                column.push_back(i + 1);
            }

            if (left != 0 && x == 0) {
                values.push_back(left);
                column.push_back(nx - 1 + nx * (y + ny * z));
            }

            if (top != 0 && y < ny - 1) {
                values.push_back(top);
                column.push_back(i + nx);
            }

            if (front != 0 && z < nz - 1) {
                values.push_back(front);
                column.push_back(i + nx * ny);
            }

            if (back != 0 && z == 0) {
                values.push_back(back);
                column.push_back(x + nx * (y + ny * (nz - 1)));
            }
        }
        row_start[nx * ny * nz] = values.size();
    }

    template<class Geometry>
	auto
    helmholtz_periodic<Geometry>::boundary(
            dt::vector<nx * nz>& top,
            dt::vector<nx * nz>& bottom) ->
		dt::vector<nx * ny * nz>
    {
        dt::vector<nx * ny * nz> b;
        for (auto i = 0u; i < nx * nz; ++i) {
            auto x = i % nx, z = i / nx;
            if (ny > 1) {
                b[x + nx * ny * z] = 8. / 3. * scale * n * n * bottom[i];
                b[x + nx * (ny - 1 + ny * z)] = 8. / 3. * scale * n * n * top[i];
            } else {
                b[x + nx * ny * z] = 3 * n * n * bottom[i];
                b[x + nx * (ny - 1 + ny * z)] = 3 * n * n * top[i];
            }
        }
        return b;
    }

    template<class Geometry>
    class helmholtz_dirichlet : public helmholtz<Geometry, Geometry::nx * (Geometry::ny-1) * Geometry::nz> {
		public:
			static constexpr unsigned int nx = Geometry::nx;
			static constexpr unsigned int ny = Geometry::ny;
			static constexpr unsigned int nz = Geometry::nz;
			static constexpr unsigned int n = Geometry::n;
        private:
            double scale;
        protected:
            using dt::matrix<nx * (ny-1) * nz>::values;
            using dt::matrix<nx * (ny-1) * nz>::column;
            using dt::matrix<nx * (ny-1) * nz>::row_start;
        public:
            dt::vector<nx * (ny-1) * nz> boundary(
                    dt::vector<nx * nz>&,
                    dt::vector<nx * nz>&);
            helmholtz_dirichlet(double);
    };

    template<class Geometry>
    helmholtz_dirichlet<Geometry>::helmholtz_dirichlet(double k)
        : helmholtz<Geometry, nx * (ny-1) * nz>(), scale(k)
    {
        for (auto i = 0u; i < nx * (ny-1) * nz; ++ i) {
            row_start[i] = values.size();
            auto x = i % nx, y = (i % (nx * (ny-1))) / nx, z = i / (nx * (ny-1));
            double diagonal = 1. - scale * 6. * (n * n),
                   offdiagonal = scale * n * n,
                   left = offdiagonal,
                   right = offdiagonal,
                   top = offdiagonal,
                   bottom = offdiagonal,
                   front = offdiagonal,
                   back = offdiagonal;

            if (nx == 1) {
                diagonal += left + right;
                left = right = 0;
            } else if (nx == 2) {
                if (x == 0) {
                    right += left;
                    left = 0;
                } else {
                    left += right;
                    right = 0;
                }
            }

            if (y == 0) { bottom = 0; }
            if (y == ny - 2) { top = 0; }

            if (nz == 1) {
                diagonal += back + front;
                back = front = 0;
            } else if (nz == 2) {
                if (z == 0) {
                    front += back;
                    back = 0;
                } else {
                    back += front;
                    front = 0;
                }
            }

            if (front != 0 && z == nz - 1) {
                values.push_back(front);
                column.push_back(x + nx * y);
            }

            if (back != 0 && z > 0) {
                values.push_back(back);
                column.push_back(i - nx * (ny - 1));
            }

            if (bottom != 0 && y > 0) {
                values.push_back(bottom);
                column.push_back(i - nx);
            }

            if (right != 0 && x == nx - 1) {
                values.push_back(right);
                column.push_back(nx * (y + (ny - 1) * z));
            }

            if (left != 0 && x > 0) {
                values.push_back(left);
                column.push_back(i - 1);
            }

            if (diagonal != 0) {
                values.push_back(diagonal);
                column.push_back(i);
            }

            if (right != 0 && x < nx - 1) {
                values.push_back(right);
                column.push_back(i + 1);
            }

            if (left != 0 && x == 0) {
                values.push_back(left);
                column.push_back(nx - 1 + nx * (y + (ny - 1) * z));
            }

            if (top != 0 && y < ny - 2) {
                values.push_back(top);
                column.push_back(i + nx);
            }

            if (front != 0 && z < nz - 1) {
                values.push_back(front);
                column.push_back(i + nx * (ny - 1));
            }

            if (back != 0 && z == 0) {
                values.push_back(back);
                column.push_back(x + nx * (y + (ny -1) * (nz - 1)));
            }
        }
        row_start[nx * (ny - 1) * nz] = values.size();
    }

    template<class Geometry>
	auto
    helmholtz_dirichlet<Geometry>::boundary(
            dt::vector<nx * nz>& top,
            dt::vector<nx * nz>& bottom) ->
		dt::vector<nx * (ny-1) * nz>
    {
        dt::vector<nx * (ny-1) * nz> b;
        if (ny > 1) {
            for (auto i = 0; i < nx * nz; ++i) {
                auto x = i % nx, z = i / nx;
                b[x + nx * (ny - 1) * z] = scale * n * n * bottom[i];
                b[x + nx * (ny - 2 + (ny - 1) * z)] = scale * n * n * top[i];
            }
        }
        return b;
    }
}

#endif
