#include "interface.h"
#include "vector.h"
#include "dense_matrix.h"
#include "sparse_matrix.h"

void
gemv(vector& y, double alpha, dense_matrix a, vector x, double beta)
{
	cblas_dgemv(CblasRowMajor, CblasNoTrans, a.rows(), a.cols(), alpha,
			a.m_values, a.cols(), x.m_values, x.cols(), beta, y.m_values, y.cols());
}

void
gemv(vector& y, double alpha, sparse_matrix a, vector x, double beta)
{
	unsigned int row = 0;
	double row_value = y[row];
	double mult_value = 0;
	for (unsigned int i = 0; i < a.num_elements(); ++i) {
		while (row < a.size() - 1 && i >= a.m_row_starts[row]) {
			y[row] = alpha * mult_value + beta * row_value;
			mult_value = 0;
			row_value = y[++row];
		}
		mult_value += a.m_values[i] * x[a.m_cols[i]];
	}
	y[row] = alpha * mult_value + beta * row_value;
}
