#ifndef TP1_METODOS_SVD_H
#define TP1_METODOS_SVD_H

#include <cassert>
#include <chrono>
#include "Matrix.h"

#define GET_TIME std::chrono::high_resolution_clock::now()
#define     GET_TIME_DELTA(begin, end) \
     std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()

tuple<Matrix, Matrix> pca(const Matrix &A, unsigned int num_components);

tuple<Matrix, double> power_method(Matrix& x_0, Matrix& input);
tuple<Matrix, double> power_method(Matrix& x_0, Matrix& input, int max_iters);

tuple<Matrix, double> powerMethodQ1(Matrix x_0, const Matrix &a);
tuple<Matrix, double> powerMethodQ1(Matrix x_0, const Matrix &a, long N);

tuple<Matrix, Matrix> svd(const Matrix &A, unsigned int num_components);

Matrix ones(int rows, int cols);

Matrix random(int rows, int cols);

#endif //TP1_METODOS_SVD_H
