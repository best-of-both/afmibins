#include <iostream>

#include "optimizations.h"
#include "interface.h"
#include "vector.h"
#include "dense_matrix.h"

vector&
operator+=(vector& v, future<future<double, dense_matrix, mult_op>, vector, mult_op> w)
{
	auto l = w.left;
	gemv(v, l.left, l.right, w.right, 1);
	return v;
}

vector&
operator+=(vector& v, future<double, future<dense_matrix, vector, mult_op>, mult_op> w)
{
	auto r = w.right;
	gemv(v, w.left, r.left, r.right, 1);
	return v;
}

vector&
operator+=(vector& v, future<dense_matrix, future<double, vector, mult_op>, mult_op> w)
{
	auto r = w.right;
	gemv(v, r.left, w.left, r.right, 1);
	return v;
}

vector&
operator+=(vector& v, future<dense_matrix, vector, mult_op> w)
{
	gemv(v, 1, w.left, w.right, 1);
	return v;
}

vector&
operator-=(vector& v, future<future<double, dense_matrix, mult_op>, vector, mult_op> w)
{
	auto l = w.left;
	gemv(v, -l.left, l.right, w.right, 1);
	return v;
}

vector&
operator-=(vector& v, future<double, future<dense_matrix, vector, mult_op>, mult_op> w)
{
	auto r = w.right;
	gemv(v, -w.left, r.left, r.right, 1);
	return v;
}

vector&
operator-=(vector& v, future<dense_matrix, future<double, vector, mult_op>, mult_op> w)
{
	auto r = w.right;
	gemv(v, -r.left, w.left, r.right, 1);
	return v;
}

vector&
operator-=(vector& v, future<dense_matrix, vector, mult_op> w)
{
	gemv(v, -1, w.left, w.right, 1);
	return v;
}
