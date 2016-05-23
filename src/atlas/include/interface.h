#ifndef ATLAS_BLAS_INTERFACE_H
#define ATLAS_BLAS_INTERFACE_H

class vector;
class dense_matrix;
class sparse_matrix;

void gemv(vector&, double, dense_matrix, vector, double);
void gemv(vector&, double, sparse_matrix, vector, double);

#endif
