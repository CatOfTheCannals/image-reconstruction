#ifndef GEOMETRIAS_H
#define GEOMETRIAS_H

#include "discretizador.h"

#include <iostream>

using namespace std;

void static aux_geometria_circular(std::pair <int,int> center, int base_index_x, int base_index_y, Discretizador discretizacion, int k) { //auxiliar para disparar una geo circular dado un punto.
	//tenemos que x^2 + y^2 = r^2
	//si usamos n^2 como r, r^2 = n^4
	unsigned tamanio_matrix = discretizacion.img_size();
	int r = round(tamanio_matrix*tamanio_matrix * k); // elegido radio  = n^2 para elegir algo arbitrariamente "grande"
	int squared_r = r*r; // puede que esté mal
	int x = -1 * r + base_index_x; // (x = sqrt(r^2 - y^2))
	int y = 0 + base_index_y;  	// (y = sqrt(r^2 - x^2))

	for (int i = 0; i < 2*r; ++i) { // ciclo x
		discretizacion.trazar_rayo(make_pair(x, y), center);//dispara el rayo
		x++;//se mueve en el eje x
		y = round(sqrt(squared_r - x^2));//se mueve en el eje y)
	}

	y = -1 * r + base_index_x; // (y = sqrt(r^2 - x^2)) //
	x = 0 + base_index_y;  	// (x = sqrt(r^2 - y^2))

	for (int i = 0; i < 2*r; ++i) { // ciclo y
		discretizacion.trazar_rayo(make_pair(x, y), center);//dispara el rayo
		y++;//se mueve en el eje y
		x = round(sqrt(squared_r - y^2));//se mueve en el eje x)
	}
}

bool static es_par(int n) {
	int x = n%2;
	return x == 0;
}

void static definir_variables_centro(int base_index, std::pair <int,int> center, int tamanio_matrix) {
	if(es_par(tamanio_matrix)) {
		center = make_pair(tamanio_matrix/2, tamanio_matrix/2);
		base_index = tamanio_matrix/2;
	} else {
		center = make_pair(tamanio_matrix/2 + 1, tamanio_matrix/2 + 1);
		base_index = tamanio_matrix/2 +1;
	}
}

//Dado una cantidad de rayos n, usa la función auxiliar geo_diagonales para trazar
//n rayos diagonales.
struct geo_diagonalesPorCantidad {
	int _n;

	//para pasar parámetros.
	geo_diagonalesPorCantidad(int n) : _n(n) {
	}

	//Aplica la geometria (genera los rayos)
	void aplicar(Discretizador &discretizacion) {
		unsigned img_size = discretizacion.img_size();
		int count = _n;
		unsigned cant_diag = 2 * img_size - 1;
		if (count > cant_diag) {
			int separacionEntreRayos = 1;
			while (count > cant_diag) {
				for (unsigned i = 1; i <= cant_diag; i += separacionEntreRayos) {
					discretizacion.trazar_rayo(make_pair(i, -1), make_pair(-1, i));
					discretizacion.trazar_rayo(make_pair(-1, img_size - 1 - i), make_pair(i, img_size));
				}
				count = count - cant_diag;
			}
		} //cuando la cantidad de rayos a tirar ya no es más grande que la cantidad de diagonales que puedo tirar,
		//tiro los rayos restantes a intervalos equidistantes
		int separacionEntreRayos = (img_size / count)*2;  //Fixme: Adri, por favor revisá como se generan los rayos a
		// partir de acá en adelante.

		for (unsigned i = 1; i <= cant_diag; i += separacionEntreRayos) {
			discretizacion.trazar_rayo(make_pair(i, -1), make_pair(-1, i));
			discretizacion.trazar_rayo(make_pair(-1, img_size - 1 - i), make_pair(i, img_size));
		}

	}
};


//Geometria que traza rayos diagonales
struct geo_diagonales{

	int _c;

	//para pasar parametros. 
	geo_diagonales(int c) : _c(c){
	}

	//Aplica la geometria (genera los rayos)
	void aplicar(Discretizador& discretizacion){
		unsigned img_size = discretizacion.img_size();
		unsigned cant_diag = 2*img_size - 1;

		for(unsigned i = 1; i <= cant_diag; i+=_c){
			discretizacion.trazar_rayo(make_pair(i, -1), make_pair(-1, i));
			discretizacion.trazar_rayo(make_pair(-1, img_size - 1 - i), make_pair(i, img_size));
		}
	}
};

//fixme: tira error que no puedo arreglar, falta adapatarla para
// fixme:que tome cantidad de rayos totales y cantidad de focos


struct geo_punto_fijo_random{

	int _cant_puntos;
	int _cant_rayos;

	//En cada borde de la imagen toma cant_puntos pixeles aleatorios y traza
	//Desde ellos cant_rayos aleatorios a cada uno de los bordes restantes (3*cant_rayos en total, no necesariamente distintos)
	geo_punto_fijo_random(int cant_puntos, int cant_rayos){
		_cant_puntos = cant_puntos;
		_cant_rayos = cant_rayos;
	}

	void aplicar(Discretizador& discretizacion){

		unsigned img_size = discretizacion.img_size();

		//Borde izq a los demas
		for(int i = 0; i < _cant_puntos; i++){
			unsigned k = rand() % img_size;
			for(int j = 0; j < _cant_rayos; j++){
				int dest = 0;
				do{ dest = rand() % img_size; } while(dest == k && k == 0);
				discretizacion.trazar_rayo(make_pair(k, 0), make_pair(0, dest)); //borde sup.
				discretizacion.trazar_rayo(make_pair(k, 0), make_pair(rand() % img_size, img_size-1)); //borde der.
				discretizacion.trazar_rayo(make_pair(k, 0), make_pair(img_size-1, rand() % img_size)); //borde inf.

			}
		}

		//Borde der a los demas
		for(int i = 0; i < _cant_puntos; i++){
			unsigned k = rand() % img_size;
			for(int j = 0; j < _cant_rayos; j++){
				int dest = 0;
				do{ dest = rand() % img_size; } while(dest == k && k == img_size-1); //solo para que no coincida el punto.
				discretizacion.trazar_rayo(make_pair(k, img_size-1), make_pair(0, dest)); //borde sup.
				discretizacion.trazar_rayo(make_pair(k, img_size-1), make_pair(rand() % img_size, 0)); //borde izq.
				discretizacion.trazar_rayo(make_pair(k, img_size-1), make_pair(img_size-1, rand() % img_size)); //borde inf.

			}
		}
	}

};

//Dado f (cantidad de los focos) y c (cantidad de rayos),
// traza los f rayos que van desde los focos de un lateral al lateral opuesto de la imagen
struct geo_focosRegularesIzquierdaDerechaFijos {

	int _c;
	int _f;

	//para pasar parámetros.
	geo_focosRegularesIzquierdaDerechaFijos(int c, int f) {
		_c= c;
		_f = f;
	}

	//Aplica la geometria (genera los rayos)
	void aplicar(Discretizador &discretizacion) {
        unsigned img_size = discretizacion.img_size();
        //f es divisor de img_size
        //c es multiplo o divisor de de img_size
        //c es multiplo de f

        int separacionFocos = img_size / ((_f/2) + (_f%2));
        for(int f = 0; f < _f; f++){

        	int rayosPorFoco = _c/_f;
        	int separacionEntreRayos = img_size / rayosPorFoco;
        	if(separacionEntreRayos != 0){
	        	for(int r = 0; r < rayosPorFoco; r++){
	        		if(f % 2 == 0)
	        			discretizacion.trazar_rayo(make_pair((f/2)*separacionFocos, 0), make_pair((r+1)*separacionEntreRayos, img_size - 1));
	        		else
	        			discretizacion.trazar_rayo(make_pair(((f/2))*separacionFocos, img_size - 1), make_pair((r+1)*separacionEntreRayos, 0));
	        	}

        	}else{
        		int countRayos = 0;
        		int actual_pixel = 0;
        		while(countRayos < rayosPorFoco){
        			if(f % 2 == 0)
	        			discretizacion.trazar_rayo(make_pair((f/2)*separacionFocos, 0), make_pair(actual_pixel, img_size - 1));
	        		else
	        			discretizacion.trazar_rayo(make_pair(((f/2))*separacionFocos, img_size - 1), make_pair(actual_pixel, 0));
        		
        			countRayos++;
        			actual_pixel++;
        			if (actual_pixel == img_size) actual_pixel = 0;
        		}
        	}
        }

    }
};


//Dado f (distancia entre los focos) y c (separación entre los extremos de los rayos),
// traza los rayos que van desde los focos de un lateral al lateral opuesto de la imagen
struct geo_focosRegularesIzquierdaDerecha {

	int _separacionRayos;
	int _distanciaFocos;

	//para pasar parámetros.
	geo_focosRegularesIzquierdaDerecha(int c, int f) {
		_separacionRayos = c;
		_distanciaFocos = f;
	}

	//Aplica la geometria (genera los rayos)
	void aplicar(Discretizador &discretizacion) {
		unsigned img_size = discretizacion.img_size();

		for (unsigned i = 0; i * _distanciaFocos < img_size; i++) {
			for (unsigned j = 0; j * _separacionRayos < img_size; j++) {

				discretizacion.trazar_rayo(make_pair(i * _distanciaFocos, 0),
										   make_pair(j * _separacionRayos, img_size - 1)); //focos desde la izquierda
				discretizacion.trazar_rayo(make_pair(i * _distanciaFocos, img_size - 1),
										   make_pair(j * _separacionRayos, 0)); //focos desde la derecha
			}
		}
	}
};



//Dado f (cantidad de los focos) y c (cantidad de rayos),
// traza los f rayos que van desde los focos desde arriba a abajo y viceversa.
//Dado f (cantidad de los focos) y c (cantidad de rayos),
// traza los f rayos que van desde los focos de un lateral al lateral opuesto de la imagen
struct  geo_focosRegularesArribaAbajoFijos {

    int _c;
    int _f;

    //para pasar parámetros.
    geo_focosRegularesArribaAbajoFijos(int c, int f) {
        _c = c;
        _f = f;
    }

    //Aplica la geometria (genera los rayos)
    void aplicar(Discretizador &discretizacion) {
        unsigned img_size = discretizacion.img_size();
        //f es divisor de img_size
        //c es multiplo o divisor de de img_size
        //c es multiplo de f

        int separacionFocos = img_size / ((_f/2) + (_f%2));
        for(int f = 0; f < _f; f++){

        	int rayosPorFoco = _c/_f;
        	int separacionEntreRayos = img_size / rayosPorFoco;
        	if(separacionEntreRayos != 0){

	        	for(int r = 0; r < rayosPorFoco; r++){
	        		if(f % 2 == 0)
	        			discretizacion.trazar_rayo(make_pair(0, (f/2)*separacionFocos), make_pair(img_size - 1, (r+1)*separacionEntreRayos));
	        		else
	        			discretizacion.trazar_rayo(make_pair(img_size - 1, ((f/2))*separacionFocos), make_pair(0, (r+1)*separacionEntreRayos));
	        	}
        	}
        	else{
        		int countRayos = 0;
        		int actual_pixel = 0;
        		while(countRayos < rayosPorFoco){
        			if(f % 2 == 0)
	        			discretizacion.trazar_rayo(make_pair(0, (f/2)*separacionFocos), make_pair(img_size - 1, actual_pixel));
	        		else
	        			discretizacion.trazar_rayo(make_pair(img_size - 1, ((f/2))*separacionFocos), make_pair(0, actual_pixel));
        		
        			countRayos++;
        			actual_pixel++;
        			if (actual_pixel == img_size) actual_pixel = 0;
        		}
        	}
        }
    }
};

//Dado f (distancia entre los focos) y c (separación entre los extremos de los rayos),
// traza los rayos que van desde los focos desde arriba de la imagen, hasta abajo y viceversa.
struct geo_focosRegularesArribaAbajo {

	int _separacionRayos;
	int _distanciaFocos;

	//para pasar parámetros.
	geo_focosRegularesArribaAbajo(int c, int f) {
		_separacionRayos = c;
		_distanciaFocos = f;
	}

	//Aplica la geometria (genera los rayos)
	void aplicar(Discretizador &discretizacion) {
		unsigned img_size = discretizacion.img_size();

		for (unsigned i = 0; i * _distanciaFocos < img_size; i++) {
			for (unsigned j = 0; j * _separacionRayos < img_size; j++) {

				discretizacion.trazar_rayo(make_pair(0, i * _distanciaFocos), make_pair(img_size - 1, j * _separacionRayos)); //focos desde arriba.
				discretizacion.trazar_rayo(make_pair(img_size - 1, i * _distanciaFocos), make_pair(0, j * _separacionRayos)); //focos desde abajo.
			}
		}
	}
};






//Dado una cantidad de rayos n, traza n rayos horizontales.
struct geo_horizontalesPorCantidad {
	int _n;

	//para pasar parámetros.
	geo_horizontalesPorCantidad(int n) : _n(n){
	}

	void aplicar(Discretizador& discretizacion) {

		//usa geo_horizontales para generar los rayos
		unsigned img_size = discretizacion.img_size();
		int count = _n;
		if (count > img_size) {
			while (count > img_size) {
				int separacionEntreRayos = 1;
				//Aplica la geometria horizontales (genera los rayos en cada píxel)
				for (unsigned i = 0; i * separacionEntreRayos < img_size; i++) {
					discretizacion.trazar_rayo(make_pair(i * separacionEntreRayos, 0),
											   make_pair(i * separacionEntreRayos, img_size - 1));
				}
				count = count - img_size;
			}
		} //cuando la cantidad de rayos a tirar ya no es más grande que la cantidad de píxeles por columna,
		//tiro los rayos restantes a intervalos equidistantes
		int separacionEntreRayos = img_size / count;
		//Aplica la geometria horizontales (genera los rayos en cada píxel)
		for (unsigned i = 0; i * separacionEntreRayos < img_size; i++) {
			discretizacion.trazar_rayo(make_pair(i * separacionEntreRayos, 0),
									   make_pair(i * separacionEntreRayos, img_size - 1));
		}
	}
};



//Dado un intervalo c (cada cuanto va a haber un rayo),
// devuelve vector de inicios y fines de rayos.
//Geometría que traza rayos horizontales separados por un intervalo c.
struct geo_horizontales {

	int _c;

	//para pasar parámetros.
	geo_horizontales(int c) : _c(c){

	}

	//Aplica la geometria (genera los rayos)
	void aplicar(Discretizador& discretizacion){
		unsigned img_size = discretizacion.img_size();

		for(unsigned i = 0; i*_c < img_size; i++){
			discretizacion.trazar_rayo(make_pair(i*_c, 0), make_pair(i*_c, img_size-1));
		}
	}



};



//Dado una cantidad de rayos n, traza n rayos verticales.
struct geo_verticalesPorCantidad {
	int _n;

	//para pasar parámetros.
	geo_verticalesPorCantidad(int n) : _n(n){
	}

	void aplicar(Discretizador& discretizacion) {

		//usa geo_horizontales para generar los rayos
		unsigned img_size = discretizacion.img_size();
		int count = _n;
		if (count > img_size) {
			while (count > img_size) {
				int separacionEntreRayos = 1;
				//Aplica la geometria verticales (genera los rayos en cada píxel)
				for (unsigned i = 0; i * separacionEntreRayos < img_size; i++) {
					discretizacion.trazar_rayo(make_pair(0, i*separacionEntreRayos),
											   make_pair(img_size-1,i*separacionEntreRayos));
				}
				count = count - img_size;
			}
		} //cuando la cantidad de rayos a tirar ya no es más grande que la cantidad de píxeles por columna,
		//tiro los rayos restantes a intervalos equidistantes
		int separacionEntreRayos = img_size / count;
		//Aplica la geometria horizontales (genera los rayos en cada píxel)
		for (unsigned i = 0; i * separacionEntreRayos < img_size; i++) {
			discretizacion.trazar_rayo(make_pair(0, i*separacionEntreRayos),
									   make_pair(img_size-1,i*separacionEntreRayos));
		}
	}
};


//Geometría que traza rayos verticales separados por un intervalo c.
struct geo_verticales {

	int _c;

	//para pasar parámetros.
	geo_verticales(int c) : _c(c){

	}

	//Aplica la geometria (genera los rayos)
	void aplicar(Discretizador& discretizacion){
		unsigned img_size = discretizacion.img_size();

		for(unsigned i = 0; i*_c < img_size; i++){
			discretizacion.trazar_rayo(make_pair(0, i*_c), make_pair(img_size-1,i*_c));
		}
	}

};



//Dado una cantidad de rayos n, traza n rayos verticales.
struct geo_centralPorCantidad {
	int _n;

	//para pasar parámetros.
	geo_centralPorCantidad(int n) : _n(n){
	}

	void aplicar(Discretizador& discretizacion) {
		//usa geo_horizontales para generar los rayos
		unsigned img_size = discretizacion.img_size();
		int count = _n;
		unsigned centro = discretizacion.img_size() / 2;
		if (count > 2 * img_size) {
			while (count > 2 * img_size) {
				int separacionEntreRayos = 1;
				//Aplica la geometria central, (desde cada píxel del borde).
				for (unsigned i = 0; i * separacionEntreRayos <
									 img_size; i++) { //rayos que atraviesan el centro y el borde superior e inferior de la imagen.
					discretizacion.trazar_rayo(make_pair(0, i * separacionEntreRayos), make_pair(centro, centro));
				}
				for (unsigned i = 0; i * separacionEntreRayos < img_size; i++) { //rayos que llegan al lateral superior
					discretizacion.trazar_rayo(make_pair(i * separacionEntreRayos, 0), make_pair(centro, centro));
				}
				count = count - 2 * img_size;
			}
		} //cuando la cantidad de rayos a tirar ya no es más grande que la cantidad de píxeles por columna*2,
		//tiro los rayos restantes a intervalos equidistantes
		int separacionEntreRayos = img_size / (count / 2);
		//Aplica la geometria central, (desde cada píxel del borde).
		for (unsigned i = 0; i * separacionEntreRayos <
							 img_size; i++) { //rayos que atraviesan el centro y el borde superior e inferior de la imagen.
			discretizacion.trazar_rayo(make_pair(0, i * separacionEntreRayos), make_pair(centro, centro));
		}
		for (unsigned i = 0; i * separacionEntreRayos < img_size; i++) { //rayos que llegan al lateral superior
			discretizacion.trazar_rayo(make_pair(i * separacionEntreRayos, 0), make_pair(centro, centro));
		}
	}
};

//Geometría que traza rayos desde el centro de la imagen hacia los cuatro laterales, separados por un intevalo c,
// versión simple de circular.
struct geo_central {

	int _c;

	//para pasar parámetros.
	geo_central (int c) : _c(c){

	}

	//Aplica la geometria (genera los rayos)
	void aplicar(Discretizador& discretizacion){
		unsigned img_size = discretizacion.img_size();
		unsigned centro = discretizacion.img_size()/2;

		for(unsigned i = 0; i*_c < img_size; i++){ //rayos que atraviesan el centro y el borde superior e inferior de la imagen.
			discretizacion.trazar_rayo(make_pair(0, i*_c), make_pair(centro,centro));
		}
		for(unsigned i = 0; i*_c < img_size; i++){ //rayos que llegan al lateral superior
			discretizacion.trazar_rayo(make_pair(i*_c, 0), make_pair(centro,centro));
		}
	}
};

//traza k rayos aleatorios
struct geo_random{

	int _k;

	geo_random(int k): _k(k){

	}

	void aplicar(Discretizador& discretizacion){
		
		unsigned img_size = discretizacion.img_size();

		for(unsigned i = 0; i < _k; i++){
			unsigned y0 = rand() % img_size;
			unsigned x0 = rand() % img_size;

			unsigned y1 = rand() % img_size;
			unsigned x1 = rand() % img_size;

			discretizacion.trazar_rayo(make_pair(y0, x0), make_pair(y1, x1));
		}
	}
};

struct geo_circular {
	double _k;//para aumentar el radio del círculo

	geo_circular(double k): _k(k){

	}

	void aplicar(Discretizador& discretizacion){

		//Matrix: 	[(0,0) (1,0) (2,0)] Tamanio: 3, n
		//			[(0,1) (1,1) (2,1)]
		//			[(0,2) (1,2) (2,2)]

		unsigned img_size = discretizacion.img_size();
		std::pair <int,int> center;
		int base_index = 0;

		definir_variables_centro(base_index, center, img_size);
		aux_geometria_circular(center, base_index, base_index, discretizacion, _k);
	}
};

struct geo_focos_circular {

	double _k;//para aumentar el radio de los círculos

	geo_focos_circular(double k): _k(k){

	}

	void aplicar(Discretizador& discretizacion) {

		unsigned img_size = discretizacion.img_size();

		int base_index_x = 0;
		int base_index_y = 0;
		//quiero que el base index esté
		//disparamos circulares en cada pixel de los bordes
		for (int i = 0; i < img_size; ++i) { //borde superior e inferior (itero x)
			aux_geometria_circular(make_pair(i, 0), 	   base_index_x, base_index_y, 			  discretizacion, _k);
			aux_geometria_circular(make_pair(i, img_size), base_index_x, base_index_y + img_size, discretizacion, _k);
			base_index_x++;
		}
		base_index_x = 0;
		for (int i = 0; i < img_size; ++i) { //borde izquierdo y derecho (itero y)
			aux_geometria_circular(make_pair(0, 	   i), base_index_x, 			base_index_y, discretizacion, _k);
			aux_geometria_circular(make_pair(img_size, i), base_index_x + img_size, base_index_y, discretizacion, _k);
			base_index_y++;
		}
	}
};

#endif

