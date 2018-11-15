#ifndef TP2_METODOS_UTILS_H
#define TP2_METODOS_UTILS_H

#include <iostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include "chrono"
#include <stdlib.h>
#include <tuple>

#define GET_TIME std::chrono::high_resolution_clock::now()
#define GET_TIME_DELTA(begin, end) \
     std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()

#define LOG_EXP_PCA_VECS false
#define LOG_KNN_DISTANCES false
#define LOG_EXP_RESULTS_SUM true
#define DATASET_DEBUG_MODE true


#endif //TP2_METODOS_UTILS_H
