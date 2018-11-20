//
// Created by juan on 16/06/18.
//

#ifndef PROJECT_RAYOS_H
#define PROJECT_RAYOS_H

#include "matrix.h"
#include <vector>


/* Devuelve secuencia de pixeles por los que pasa un rayo, tener en cuenta que un pixel se considera de altura 1 */
void pixelesQueAtraviezaR(const Matrix& M, double A, double b);


#endif //PROJECT_RAYOS_H
