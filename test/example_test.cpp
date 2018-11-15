#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "gtest/gtest.h"

// --------- SET UP --------------
class runTest : public ::testing::Test {
protected:
    virtual void SetUp() {
        int i = 0;
        while (i < 10){
            int j = 0;
            while (j < 10) {
                bigArray[i][j] = rand() % 100;//Array con elementos random
                j++;
            }
            i++;
        }
        i = 0;
        while (i < 1000){
            reallyBig[i] = rand() % 200;
            i++;
        }
    }
    int bigArray[10][10]; //100 arrays de 10 elementos random
    int reallyBig[1000];
    int naturals[10];
};

TEST_F (runTest, nombreTest){


    int i = 0;
    while (i < 10) {
        naturals [i] = i;
        i++;
    }



}