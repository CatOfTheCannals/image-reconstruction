//
// Created by juan on 16/06/18.
//
#ifndef MATRIX_ALGOS_CPP
#define MATRIX_ALGOS_CPP

#include "matrix_algos.h"

/* Devuelv secuencia de pixeles por los que pasa un rayo */
void pixelesQueAtraviezaR(const Matrix& M, double A, double b) {
    std::vector <std::tuple<int, int>> res;
    std::tuple<int,double> aux;
    for (int x = 0; x < (int)(M.cols()); ++x) { //para cada bodre de cada columna calcúlo "y".
        double ydecimal;
        double yAnteriorDecimal;
        if (res.size() > 0){
            yAnteriorDecimal = ydecimal; //me guardo el valor decimal del último "y"
        }
        double y = A*x + b;
        ydecimal = y; //en la siguiente línea voy a truncar "y", por eso me guardo la versión double.
        y = (int)(y);
        if (y <= M.rows()){  //chequeo que "y" está dentro de la imagen
            if (y != ydecimal) { // chequeo que no pase por el vértice del píxel.
                res.push_back(std::make_tuple(x,y));
                res.push_back(std::make_tuple(x + 1,y));
                aux = std::make_tuple(x ,ydecimal);
            } else { // "y" pasa exactamente el vértice (parte entera de "y" = "y")
                if (res.size() > 0) {
                    int ultimoYDeRes = std::get<1>(res[res.size()-1]);
                    if(ultimoYDeRes == yAnteriorDecimal) { //solamente me importa cuando el anterior pasa por la diagonal.
                        if (ultimoYDeRes < y) { //es una diagonal ascendente.
                            res.push_back(std::make_tuple(x - 1, y - 1));
                        } else if (ultimoYDeRes > y) { //es una diagonal descendente.
                            res.push_back(std::make_tuple(x - 1, y));
                        } else { // ultimoYDeRes = y, si pasan por la arista de los píxeles, pinto el pixel de abajo.
                            if (y == 0) {
                                res.push_back(std::make_tuple(x - 1, y));
                            } else {
                                res.push_back(std::make_tuple(x - 1, y - 1));
                            }
                        }
                    }
                }
                aux = std::make_tuple(x ,ydecimal);//inclusosi no tenía valores previos guardados en el vector, me guardo aux
                // para comparar en el siguiente paso.
            }
        }
    }
}


#endif
