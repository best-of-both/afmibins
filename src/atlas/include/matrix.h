#ifndef ATLAS_MATRIX_H
#define ATLAS_MATRIX_H

#include <cstddef>

typedef struct matrix_size {
	const std::size_t rows, cols;
	matrix_size(std::size_t rows, std::size_t cols)
		: rows(rows), cols(cols) {}
	bool operator==(matrix_size);
	matrix_size operator*(matrix_size);
	matrix_size operator+(matrix_size);
	matrix_size operator-(matrix_size);
	int operator+(int n) { return rows * cols + n; }
	int operator-(int n) { return rows * cols - n; }
} matrix_size;

bool operator<(std::size_t, matrix_size);

class base_matrix {
	public:
		typedef std::size_t index_type;
		typedef std::size_t size_type;
	private:
		index_type *m_refs;
	protected:
		matrix_size m_size;
		double *m_values;

		bool is_last_ref() { return *m_refs == 1; }
	public:
		size_type rows() { return m_size.rows; }
		size_type cols() { return m_size.cols; }
		matrix_size size() { return m_size; }

		base_matrix(size_type, size_type);
		base_matrix(size_type, size_type, double*);
		base_matrix(const base_matrix&);
		~base_matrix();

	friend class vector;
};

#endif
