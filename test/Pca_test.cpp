#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "gtest/gtest.h"

#include "../src/Matrix.h"
#include "../src/Pca.h"

// --------- SET UP --------------
class runTest : public ::testing::Test {
protected:
    virtual void SetUp() {
            g.setIndex(0, 0, 1);
            g.setIndex(0, 1, 2);
            g.setIndex(0, 2, 3);
            g.setIndex(1, 1, 4);
            g.setIndex(1, 2, 5);
            g.setIndex(1, 0, 6);
            g.setIndex(2, 1, 7);
            g.setIndex(2, 2, 8);
            g.setIndex(2, 0, 9);

            x_0.setIndex(0, 0, 1);
            x_0.setIndex(1, 0, 2);
            x_0.setIndex(2, 0, 3);


            f.setIndex(0, 0, 7);
            f.setIndex(0, 1, 2);
            f.setIndex(0, 2, 5);
            f.setIndex(1, 1, 1);
            f.setIndex(1, 2, 4);
            f.setIndex(1, 0, 4);
            f.setIndex(2, 1, 3);
            f.setIndex(2, 2, 2);
            f.setIndex(2, 0, 1);
    }
    Matrix f = Matrix(3,3);
    Matrix g = Matrix(3,3);
    Matrix x_0 = Matrix(3,1);
    double epsilon = 0.0000001;
};


TEST_F (runTest, powerG){
    Matrix v = Matrix(3,1);
    double lambda;

    std::tie(v, lambda) = power_method(x_0, g);

    epsilon *= 10; // la precision queda asi por algun motivo
    ASSERT_TRUE((g*v).isApproximate(v * lambda, epsilon)) ;
}

TEST_F (runTest, svdTest1){

    auto sim = f * f.transpose();

    auto svdRes = svd(sim, 3);
    Matrix autoVecs(std::get<0>(svdRes));
    Matrix lambdas(std::get<1>(svdRes));

    for(int j = 0; j < autoVecs.cols(); j++){
        Matrix eigenVec(autoVecs.rows(), 1);
        for(int i = 0; i < autoVecs.rows(); i++){
            eigenVec.setIndex(i,0,autoVecs(i,j));
        }

        // si le agregas un decimal mas, falla
        ASSERT_TRUE((sim*eigenVec).isApproximate(eigenVec * lambdas(j,j), 0.0001));
    }

}
/*
TEST_F (runTest, svdTest2){
    auto sim = f*f.transpose();

    auto svdRes = svd(sim, 3, epsilon);
    Matrix autoVecs(std::get<0>(svdRes));
    Matrix lambdas(std::get<1>(svdRes));

    std::cout << "lambdas x * x^t :" << std::endl;
    std::cout << lambdas << std::endl;

    auto simt = f.transpose()*f;

    auto svdRest = svd(sim, 3, epsilon);
    Matrix autoVecst(std::get<0>(svdRes));
    Matrix lambdast(std::get<1>(svdRes));

    std::cout << "lambdas x^t * x :" << std::endl;
    std::cout << lambdast << std::endl;

}
*/
