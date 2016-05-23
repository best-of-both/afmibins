#include "sparse_matrix.h"
#include "futures.h"

sparse_matrix::sparse_matrix(size_type rows, size_type cols)
	: base_matrix(rows, cols) {}

future<sparse_matrix, vector, mult_op>
mult_op<sparse_matrix, vector>::wrap(sparse_matrix A, vector x)
{
	return future<sparse_matrix, vector, mult_op>(A, x);
}

sparse_matrix::~sparse_matrix()
{
	if (is_last_ref()) {
		delete[] m_row_starts;
		delete[] m_cols;
	}
}

void
mult_op<sparse_matrix, vector>::apply(result_type& y, sparse_matrix A, vector x)
{
	gemv(y, 1, A, x, 0);
}

future<double, sparse_matrix, mult_op>
mult_op<double, sparse_matrix>::wrap(double alpha, sparse_matrix x)
{
	return future<double, sparse_matrix, mult_op>(alpha, x);
}

future<sparse_matrix, sparse_matrix, add_op>
add_op<sparse_matrix, sparse_matrix>::wrap(sparse_matrix left, sparse_matrix right)
{
	assert(left.cols() == right.cols() && left.rows() == right.rows());
	return future<sparse_matrix, sparse_matrix, add_op>(left, right);
}

void
add_op<sparse_matrix, sparse_matrix>::apply(result_type& output, sparse_matrix left, sparse_matrix right)
{
	for (std::size_t i = 0; i < right.rows() * right.cols(); ++i)
		output.m_values[i] = left.m_values[i] + right.m_values[i];
}

future<sparse_matrix, sparse_matrix, sub_op>
sub_op<sparse_matrix, sparse_matrix>::wrap(sparse_matrix left, sparse_matrix right)
{
	assert(left.cols() == right.cols() && left.rows() == right.rows());
	return future<sparse_matrix, sparse_matrix, sub_op>(left, right);
}

void
sub_op<sparse_matrix, sparse_matrix>::apply(result_type& output, sparse_matrix left, sparse_matrix right)
{
	for (std::size_t i = 0; i < right.rows() * right.cols(); ++i)
		output.m_values[i] = left.m_values[i] - right.m_values[i];
}

void
mult_op<double, sparse_matrix>::apply(result_type& y, double alpha, sparse_matrix x)
{
	for (std::size_t i = 0; i < x.rows() * x.cols(); ++i)
		y.m_values[i] = alpha * x.m_values[i];
}

std::ostream&
operator<<(std::ostream& out, sparse_matrix&)
{
	out << "[]";
	return out;
}
