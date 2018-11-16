#include "Svd.h"


tuple<Matrix, Matrix> pca(const Matrix &A, unsigned int num_components) {
    Matrix X(A);
    Matrix mean(1, A.cols());

    // get "mean row"
    for(int i = 0; i < A.rows(); i++){
        mean = mean + (X.getRow(i) / X.rows());
    }
    // subtract mean row to every row: center them
    for(int i = 0; i < X.rows(); i++){
        auto centered_row = X.getRow(i) - mean;
        centered_row = centered_row / sqrt(X.rows() - 1);
        for (int j = 0; j < X.cols(); ++j) {
            X.setIndex(i,j, centered_row(j));
        }
    }

    //Matrix trans = X.transpose();
    auto M = X.mt_times_m();


    return svd(M, num_components);

}

tuple<Matrix, Matrix> svd(const Matrix &A, unsigned int num_components) {

    Matrix X(A);

    assert((num_components <= X.rows()) && (num_components <= X.cols()));

    Matrix lambdas(num_components, num_components);
    Matrix k_eigen_vectors(X.rows(), num_components);

    for (auto i = 0; i < num_components; i++) {

        Matrix x_0( random(A.rows(), 1));
        Matrix eigen_vector(A.rows(), 1);
        double eigen_value;
        // compute i_th eigen vector and its value
        tie(eigen_vector, eigen_value) = powerMethodQ1(x_0, X);
        cout << i << " eigen value: " << eigen_value << endl;

        if ((i > 0) && eigen_vector.isApproximate(k_eigen_vectors.subMatrix(0, A.rows()-1, i-1, i-1), 0.01)){
            std::cout << "se esta repitiendo autovector ñeri" << std::endl;
        }

        for (auto q = 0; q < k_eigen_vectors.rows(); q++) {
            k_eigen_vectors.setIndex(q, i, eigen_vector(q, 0)); // fill eigen vector in res matrix
        }
        lambdas.setIndex(i, i, eigen_value);
        if(i%100==0){ cout << "PCA i : " << i << endl;}
        auto external = eigen_vector*eigen_vector.transpose();
        X = X - (external * eigen_value);
    }

    return make_tuple(k_eigen_vectors, lambdas);
}

tuple<Matrix, double> power_method(Matrix& x_0, Matrix& input) {
    return power_method(x_0, input, 100);

};

tuple<Matrix, double> powerMethodQ1(Matrix x_0, const Matrix &a) {
    return powerMethodQ1(x_0, a, 25);
};

tuple<Matrix, double> powerMethodQ1(Matrix x_0, const Matrix &a, long N) {


    for (long i = 0; i < N; i++) {
        x_0 = a * x_0;
        x_0 = x_0 / x_0.norm();
    }

    Matrix producto = a * x_0;
    double val = 0.0;
    for(int i = 0 ; i< x_0.size() ; i++){
        if(x_0(i)!=0.0){
            val = producto(i)/x_0(i);
        }

    }
    std::tuple<Matrix, double> res = std::make_tuple(x_0, val);
    return res;
}

tuple<Matrix, double> power_method(Matrix& x_0, Matrix& input, int max_iters) {
    assert(input.rows() == input.cols());
    assert(input.rows() == x_0.rows());
    assert(1 == x_0.cols());

    Matrix x(x_0 / x_0.norm());
    Matrix A(input);

    double lambda;
    double prev_lambda = 0;

    Matrix Ax = A*x;

    for(int i = 0; i < max_iters; i++) {
        double Ax_norm = Ax.norm();
        if(i==0){
            cout <<"ax norm: "<< Ax_norm <<endl; }
        x = Ax / Ax_norm;
        Ax = A*x;

        lambda = (x.transpose()*Ax)(0);
        prev_lambda = lambda;
    }
    return make_tuple(x, lambda);
}

Matrix ones(int rows, int cols)  {
    Matrix res(rows, cols);
    for (std::size_t i = 0; i < rows; i++) {
        for (std::size_t j = 0; j < cols; j++) {
            res.setIndex(i, j, 1);
        }
    }
    return res;
}

Matrix random(int rows, int cols) {
    Matrix res(rows, cols);
    srand (time(NULL));
    for (std::size_t i = 0; i < rows; i++) {
        for (std::size_t j = 0; j < cols; j++) {
            res.setIndex(i, j, rand() % 100);
        }
    }
    return res;
}
