#include <iostream>

#include "vector.h"
#include "sparse_matrix.h"
#include "optimizations.h"

int
main(void)
{
	const std::size_t N = 100; //100000;
	sparse_matrix A(N, N);
	vector x(N, 1);
	vector y(N);
	y -= A * (x / (double) N);
	std::cout << y << std::endl;
	return 0;
}
