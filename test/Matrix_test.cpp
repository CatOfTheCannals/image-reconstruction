#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "gtest/gtest.h"

#include "../src/Matrix.h"

// --------- SET UP --------------
class mockMatrices : public ::testing::Test {
protected:
    virtual void SetUp() {

        e.setIndex(0, 1, 4);
        e.setIndex(0, 2, 5);
        e.setIndex(0, 0, 6);
        e.setIndex(1, 1, 1);
        e.setIndex(1, 2, 2);
        e.setIndex(1, 0, 3);
        e.setIndex(2, 1, 4);
        e.setIndex(2, 2, 5);
        e.setIndex(2, 0, 6);
        e.setIndex(3, 1, 7);
        e.setIndex(3, 2, 8);
        e.setIndex(3, 0, 9);
        e.setIndex(4, 1, 1);
        e.setIndex(4, 2, 2);
        e.setIndex(4, 0, 3);

        m.setIndex(1, 1, 1);
        m.setIndex(1, 2, 2);
        m.setIndex(1, 0, 0);
        m.setIndex(2, 1, 4);
        m.setIndex(2, 2, 5);
        m.setIndex(2, 0, 6);
        m.setIndex(0, 1, 7);
        m.setIndex(0, 2, 8);
        m.setIndex(0, 0, 9);

        l.setIndex(1, 1, 1);
        l.setIndex(1, 2, 1);
        l.setIndex(1, 0, 1);
        l.setIndex(2, 1, 1);
        l.setIndex(2, 2, 1);
        l.setIndex(2, 0, 1);
        l.setIndex(0, 1, -8);
        l.setIndex(0, 2, 1);
        l.setIndex(0, 0, 1);

        g.setIndex(0, 0, 1);
        g.setIndex(0, 1, 2);
        g.setIndex(0, 2, 3);
        g.setIndex(1, 1, 4);
        g.setIndex(1, 2, 5);
        g.setIndex(1, 0, 6);
        g.setIndex(2, 1, 7);
        g.setIndex(2, 2, 8);
        g.setIndex(2, 0, 9);

        a.setIndex(0, 0, 1);
        a.setIndex(1, 0, 2);
        a.setIndex(2, 0, 3);

        b.setIndex(1, 0, 4);
        b.setIndex(1, 0, 5);
        b.setIndex(2, 0, 6);

        h.setIndex(0, 0, 1);
        h.setIndex(1, 0, 2);
        h.setIndex(2, 0, 3);
        h.setIndex(1, 1, 4);
        h.setIndex(1, 1, 5);
        h.setIndex(2, 1, 6);


        std::cout << "Test Matrix initialized:" << std::endl;
        std::cout << "e:" << std::endl;
        std::cout << e << std::endl;
        std::cout << "m:" << std::endl;
        std::cout << m << std::endl;
        std::cout << "l:" << std::endl;
        std::cout << l << std::endl;

    }
    Matrix a = Matrix(3,1);
    Matrix b = Matrix(3,1);
    Matrix h = Matrix(3,2);
    Matrix e = Matrix(5,3);
    Matrix m = Matrix(3,3);
    Matrix l = Matrix(3,3);
    Matrix g = Matrix(3,3);
};

TEST_F (mockMatrices, getIndex){
    ASSERT_EQ(1, e(1,1));
}

TEST_F (mockMatrices, copyConstructor){

    Matrix e_copy(e);
    for (std::size_t i = 0; i < e.rows(); i++) {
        for (std::size_t j = 0; j < e.cols(); j++) {
            ASSERT_EQ(e_copy(i, j), e(i, j));
        }
    }

}

TEST_F (mockMatrices, sumAndScalarOperator){
    auto e_copy = e;
    ASSERT_EQ(e * 2,  e + e_copy);
}

TEST_F (mockMatrices, getRow){
    int rows, cols;
    std::tie(rows, cols) = e.shape();
    for (std::size_t i = 0; i < e.rows(); i++) {
        auto e_row = e.getRow(i);
        for (std::size_t j = 0; j < e.cols(); j++) {
            ASSERT_EQ(e_row(j),  e(i, j));
        }
    }
}

TEST_F (mockMatrices, shape){
    int rows, cols;
    std::tie(rows, cols) = e.shape();
    ASSERT_EQ(rows, e.rows());
    ASSERT_EQ(cols, e.cols());

}

TEST_F (mockMatrices, equality){
    auto e_copy = e;
    ASSERT_EQ(e_copy , e);
}

TEST_F (mockMatrices, swapRows){
    int i1 = 0;
    int i2 = 4;
    auto e_row_1 = e.getRow(i1);
    auto e_row_2 = e.getRow(i2);
    e.swapRows(i1, i2);
    ASSERT_EQ(e_row_1, e.getRow(i2));
    ASSERT_EQ(e_row_2, e.getRow(i1));
}

TEST_F (mockMatrices, subMatrixRows){ // FIXME: este test solo checkea que la funcion subMatrix obtiene rows adecuadamente
    int rows, cols;
    std::tie(rows, cols) = e.shape();
    for (std::size_t i = 0; i < e.rows(); i++) {
        auto e_row = e.getRow(i);
        auto e_sub = e.subMatrix(i,i,0, cols-1);
        ASSERT_EQ(e_row, e_sub);
    }
}

TEST_F (mockMatrices, matrixMultipOperator){
    Matrix r = e*m;
    ASSERT_EQ(r.rows(),e.rows()); //TODO: check if the test ends with this assertion
    ASSERT_EQ(r.cols(),m.cols());
    //ASSERT_EQ(9*6+5+6,r(0,0)); //TODO:check if this possition is calculated correctly.
    ASSERT_EQ(3*7+1+8,r(1,1));

    int rows, cols;
    std::tie(rows, cols) = m.shape();
    Matrix id = Matrix::identity(rows);
    ASSERT_EQ(m, m*id);
}

TEST_F (mockMatrices, maxCoeff){
    int x = std::get<0>(e.maxCoeff());
    int y = std::get<1>(e.maxCoeff());
    ASSERT_EQ(3, x);
    ASSERT_EQ(0, y);
}


TEST_F (mockMatrices, abs){
    Matrix l_copy(l);
    l_copy.setIndex(0, 1, 8);
    ASSERT_EQ(l_copy, l.abs());
}

TEST_F (mockMatrices, mt_times_m){

    ASSERT_TRUE(e.mt_times_m().isApproximate(e.transpose()*e, 0.01));
    ASSERT_TRUE(e.transpose().mt_times_m().isApproximate(e*e.transpose(), 0.01));
}

TEST_F (mockMatrices, hstack){
    ASSERT_EQ(h, Matrix::hstack(a, b));
}

TEST_F (mockMatrices, saveMatrixToCsv){
    string filename("../../test/debugging_matrices/matrix");
    Matrix::saveMatrixToCsv(e, filename);
    ASSERT_EQ(Matrix::loadMatrixFromFile(filename), e);

}
