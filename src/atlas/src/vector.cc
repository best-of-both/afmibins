#include "vector.h"

vector::vector(size_type rows, double v) :
	base_matrix(rows, 1, new double[rows])
{
	for (index_type i = 0; i < rows; ++i)
		m_values[i] = v;
}

vector&
vector::operator+=(vector other)
{
	assert(size() == other.size());
	for (std::size_t i = 0; i < size(); ++i)
		operator[](i) += other[i];
	return *this;
}

vector&
vector::operator-=(vector other)
{
	assert(size() == other.size());
	for (std::size_t i = 0; i < size(); ++i)
		operator[](i) -= other[i];
	return *this;
}

future<double, vector, mult_op>
mult_op<double, vector>::wrap(double alpha, vector x)
{
	return future<double, vector, mult_op>(alpha, x);
}

void
mult_op<double, vector>::apply(result_type& y, double alpha, vector x)
{
	for (std::size_t i = 0; i < x.rows(); ++i)
		y[i] = alpha * x[i];
}

future<vector, vector, add_op>
add_op<vector, vector>::wrap(vector left, vector right)
{
	assert(left.rows() == right.rows());
	return future<vector, vector, add_op>(left, right);
}

void
add_op<vector, vector>::apply(result_type& output, vector left, vector right)
{
	for (std::size_t i = 0; i < right.rows(); ++i)
		output[i] = left[i] + right[i];
}

future<vector, vector, sub_op>
sub_op<vector, vector>::wrap(vector left, vector right)
{
	assert(left.rows() == right.rows());
	return future<vector, vector, sub_op>(left, right);
}

void
sub_op<vector, vector>::apply(result_type& output, vector left, vector right)
{
	for (std::size_t i = 0; i < right.rows(); ++i)
		output[i] = left[i] - right[i];
}

std::ostream&
operator<<(std::ostream& out, vector& v)
{
	out << "[";
	for (std::size_t i = 0; i < v.rows(); ++i) {
		out << v[i];
		if (i < v.rows() - 1)
			out << " ";
	}
	out << "]'";
	return out;
}
