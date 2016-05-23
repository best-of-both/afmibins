#include "dense_matrix.h"
#include "futures.h"

dense_matrix::dense_matrix(size_type rows, size_type cols, double v)
	: base_matrix(rows, cols, new double[rows * cols])
{
	assert(m_values != nullptr);
	for (index_type i = 0; i < size(); ++i)
		m_values[i] = v;
}

future<dense_matrix, vector, mult_op>
mult_op<dense_matrix, vector>::wrap(dense_matrix A, vector x)
{
	return future<dense_matrix, vector, mult_op>(A, x);
}

void
mult_op<dense_matrix, vector>::apply(result_type& y, dense_matrix A, vector x)
{
	gemv(y, 1, A, x, 0);
}

future<double, dense_matrix, mult_op>
mult_op<double, dense_matrix>::wrap(double alpha, dense_matrix x)
{
	return future<double, dense_matrix, mult_op>(alpha, x);
}

future<dense_matrix, dense_matrix, add_op>
add_op<dense_matrix, dense_matrix>::wrap(dense_matrix left, dense_matrix right)
{
	assert(left.cols() == right.cols() && left.rows() == right.rows());
	return future<dense_matrix, dense_matrix, add_op>(left, right);
}

void
add_op<dense_matrix, dense_matrix>::apply(result_type& output, dense_matrix left, dense_matrix right)
{
	for (std::size_t i = 0; i < right.rows() * right.cols(); ++i)
		output.m_values[i] = left.m_values[i] + right.m_values[i];
}

future<dense_matrix, dense_matrix, sub_op>
sub_op<dense_matrix, dense_matrix>::wrap(dense_matrix left, dense_matrix right)
{
	assert(left.cols() == right.cols() && left.rows() == right.rows());
	return future<dense_matrix, dense_matrix, sub_op>(left, right);
}

void
sub_op<dense_matrix, dense_matrix>::apply(result_type& output, dense_matrix left, dense_matrix right)
{
	for (std::size_t i = 0; i < right.rows() * right.cols(); ++i)
		output.m_values[i] = left.m_values[i] - right.m_values[i];
}

void
mult_op<double, dense_matrix>::apply(result_type& y, double alpha, dense_matrix x)
{
	for (std::size_t i = 0; i < x.rows() * x.cols(); ++i)
		y.m_values[i] = alpha * x.m_values[i];
}

std::ostream&
operator<<(std::ostream& out, dense_matrix&)
{
	out << "[]";
	return out;
}
