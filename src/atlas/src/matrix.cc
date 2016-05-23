#include <cassert>
#include "matrix.h"

bool
matrix_size::operator==(matrix_size other) {
	return other.rows == rows &&
		   other.cols == cols;
}

matrix_size
matrix_size::operator*(matrix_size other) {
	assert(other.rows == cols);
	return {rows, other.cols};
}

matrix_size
matrix_size::operator+(matrix_size other) {
	assert(*this == other);
	return {rows, cols};
}

matrix_size
matrix_size::operator-(matrix_size other) {
	assert(*this == other);
	return {rows, cols};
}

bool
operator<(std::size_t other, matrix_size size)
{
	return other < size.rows * size.cols;
}

base_matrix::base_matrix(size_type rows, size_type cols) :
	m_refs(new index_type(1)), m_size(rows, cols), m_values(nullptr) {}

base_matrix::base_matrix(size_type rows, size_type cols, double* values) :
	m_refs(new index_type(1)), m_size(rows, cols), m_values(values) {}

base_matrix::base_matrix(const base_matrix& other) :
	m_refs(other.m_refs), m_size(other.m_size), m_values(other.m_values)
{
	++*m_refs;
}

base_matrix::~base_matrix()
{
	if (is_last_ref()) {
		delete[] m_values;
		delete m_refs;
	} else {
		--*m_refs;
	}
}
