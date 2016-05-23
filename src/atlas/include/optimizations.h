#ifndef ATLAS_BLAS_OPTIMIZATIONS_H
#define ATLAS_BLAS_OPTIMIZATIONS_H

#include "futures.h"
#include "vector.h"
#include "dense_matrix.h"

vector& operator+=(vector&, future<future<double, dense_matrix, mult_op>, vector, mult_op>);
vector& operator+=(vector&, future<double, future<dense_matrix, vector, mult_op>, mult_op>);
vector& operator+=(vector&, future<dense_matrix, future<double, vector, mult_op>, mult_op>);
vector& operator+=(vector&, future<dense_matrix, vector, mult_op>);
vector& operator-=(vector&, future<future<double, dense_matrix, mult_op>, vector, mult_op>);
vector& operator-=(vector&, future<double, future<dense_matrix, vector, mult_op>, mult_op>);
vector& operator-=(vector&, future<dense_matrix, future<double, vector, mult_op>, mult_op>);
vector& operator-=(vector&, future<dense_matrix, vector, mult_op>);

#endif
