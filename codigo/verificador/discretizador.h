	
#ifndef DISCRETIZADOR_H
#define DISCRETIZADOR_H

#include "matrix.h"
#include "matrix_algos.h"

#include <string>
#include <stdexcept>
#include <math.h>
#include <chrono>

class Discretizador{

private:

	Matrix _img;
	std::vector< std::pair< std::pair<int, int>, std::pair<int, int> > > _rayos;
	size_t _n;
	double _ruido;
	double _tiempo_matrices;
	double _tiempo_cml;

public:

	/*
	Discretiza una imagen de filename. La imagen debe ser cuadrada y n debe ser un multiplo del ancho/altura de la img.
	*/
	Discretizador(const std::string& filename, size_t n, double r){

		_n = n;
		_img = load_matrix(filename);
		_ruido = r;

		if(_img.rows() != _img.cols()){
			_img = get_zero_initialized(n);
			throw std::runtime_error("La imagen no es cuadrada");
		}
		
		if( (_img.rows() % _n) != 0 ){
			_n = 1;
			_img = get_zero_initialized(1);
			throw std::runtime_error("error: Discretizacion invalida");
		}

	}

	void trazar_rayo(int init_i, int init_j, int end_i, int end_j){

		auto init = std::make_pair(init_i, init_j);
		auto end  = std::make_pair(end_i, end_j);
		trazar_rayo(init, end); 
	}

	void trazar_rayo(const std::pair<int, int>& init, const std::pair<int, int>& end){

		auto rayo = std::make_pair(init, end);
		_rayos.push_back(rayo);
	}

	size_t img_size(){
		return _img.rows();
	}

	Matrix ejecutar_analisis_tomografico(){

		std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();
		auto matrices = matrices_del_sistema(); //CALCULAR D y t
	  	std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
	  	_tiempo_matrices = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

		Matrix& D = matrices.first;
		Matrix& t = matrices.second;

		agregar_ruido(t);
		Matrix Dt = trasponer(D);

		begin = std::chrono::system_clock::now();
		Matrix x = resolver_sistema(Dt*D, Dt*t); //EJECUTAR CML
	  	end = std::chrono::system_clock::now();
		_tiempo_cml = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
	  	
	  	return x;

	}

	void verificar_trazado(std::string output_img, bool all_debug){


		int n = _img.rows();
		if(_rayos.size() == 0) throw std::runtime_error("No se han trazado rayos");
		

		//Cambiar img a RGB
		char* data = new char[_img.rows()*_img.cols()*3];
		char* data_ptr = data;
		for(size_t i = 0; i < _img.rows(); i++){
			for(size_t j = 0; j < _img.rows(); j++){
				data_ptr[0] = _img(i,j);
				data_ptr[1] = _img(i,j);
				data_ptr[2] = _img(i,j);
				data_ptr += 3;
			}
		}

		std::cout << "debug: _rayos.size() " << _rayos.size() << std::endl;
		//Trazar Rayos!

		auto resultado = std::make_pair(get_zero_initialized(_rayos.size(), _n*_n), get_initialized(_rayos.size(), 1, double(0)));

		//trazar rayo:
		//y = redondear(((x-x0)/(x1-x0))*(y1-y0) + y0) (recta que pasa por 2 puntos (x0, y0), (x1, y1) con su img en el conjunto N)
		for(size_t i = 0; i < _rayos.size(); i++){

			auto rayo = _rayos[i];
			auto init = rayo.first;
			auto end  = rayo.second;

			int y0 = init.first;
			int x0 = init.second;
			int y1 = end.first;
			int x1 = end.second;
			
			double x_diff = double(x1) - x0; 
			double y_diff = double(y1) - y0;

			if(x_diff == 0){ //recta vertical (pendiente infinita, es decir x_diff == 0)
				for(size_t y = 0; y < _img.rows(); y++){
					//x0 == x1

					data[(y*n + x0)*3 + 0] = 255;
					data[(y*n + x0)*3 + 1] = 0;
					data[(y*n + x0)*3 + 2] = 0;
					
				}
			}
			else{

				for(size_t x = 0; x < _img.cols(); x++){ // _img es cuadrada, es lo mismo rows y cols
					
					int y = (int) round( ((double(x)-x0)/x_diff)*y_diff + y0 );
					if (y < 0  || y >= (long int)_img.rows()){
						if(all_debug){
							std::cout << "Salteado: " << y << std::endl;
							std::cout << "\t original: " << round( ((double(x)-x0)/x_diff)*y_diff + y0 ) << std::endl;
							std::cout << "\t x:      " << x << std::endl;
							std::cout << "\t x0:     " << x0 << std::endl;
							std::cout << "\t y0:     " << y0 << std::endl;
							std::cout << "\t x_diff: " << x_diff << std::endl;
							std::cout << "\t y_diff: " << y_diff << std::endl;
							std::cout << "\t k = x-x0: " << x-x0 << std::endl;
							std::cout << "\t k =/x_diff: " << (double(x)-x0)/x_diff << std::endl;
							std::cout << "\t k =* y_diff: " << ((double(x)-x0)/x_diff)*y_diff << std::endl;
							std::cout << "\t k += y0: " << ((double(x)-x0)/x_diff)*y_diff + y0 << std::endl;
							std::cout << "\t round(k): " << round( ((double(x)-x0)/x_diff)*y_diff + y0 ) << std::endl;
						}
						continue; //No sirve ningún pixel por fuera de la discretización.
					}
			
					data[(y*n + x)*3 + 0] = 255;
					data[(y*n + x)*3 + 1] = 0;
					data[(y*n + x)*3 + 2] = 0;
				}
			}
			
		}

		if(!SavePPMFile(output_img.c_str(), data, n, n, PPM_LOADER_PIXEL_TYPE_RGB_8B, NULL)){
			throw std::runtime_error("Error guardando imagen");	
		} 
		delete[] data;

		return;
	}

	double tiempo_armado_matrices(){
		return _tiempo_matrices;
	}

	double tiempo_cml(){
		return _tiempo_cml;
	}

private:

	//Agrega ruido a t. 
	//requiere: t es un vector columna
	void agregar_ruido(Matrix& t){

		srand(0); //resultados reproducibles;

		double fMax = _ruido;
		double fMin = 0.0; 

		for(size_t i = 0; i < t.rows(); i++){
			double f = (double)rand() / RAND_MAX;
    		double rand_double = fMin + f * (fMax - fMin);
			t(i, 0) = t(i,0) + rand_double;
		}
	}

	//Devuelve un par <D, t>, donde D es la matriz de distancias y t el vector de velocidades, ambos construidos con los rayos trazados.
	//D tiene dimension m x n^2, donde m es la cantidad de rayos y n el tamaño de la discretizacion
	//t tiene dimension m x 1, donde m es la cantidad de rayos.
	std::pair<Matrix, Matrix> matrices_del_sistema(){

		if(_rayos.size() == 0) throw std::runtime_error("No se han trazado rayos");
		
		auto resultado = std::make_pair(get_zero_initialized(_rayos.size(), _n*_n), get_initialized(_rayos.size(), 1, double(0)));

		Matrix& D = resultado.first;
		Matrix& t = resultado.second;
		std::vector<Matrix::value_type> velocidades_rayo_actual(_n*_n, Matrix::value_type(0));
		std::vector<Matrix::value_type> num_vel(_n*_n, Matrix::value_type(0));

		//trazar rayo:
		//y = redondear(((x-x0)/(x1-x0))*(y1-y0) + y0) (recta que pasa por 2 puntos (x0, y0), (x1, y1) con su img en el conjunto N)
		for(size_t i = 0; i < _rayos.size(); i++){

			auto rayo = _rayos[i];
			auto init = rayo.first;
			auto end  = rayo.second;

			int y0 = init.first;
			int x0 = init.second;
			int y1 = end.first;
			int x1 = end.second;
			
			double x_diff = double(x1) - x0; 
			double y_diff = double(y1) - y0;

			if(x_diff == 0){ //recta vertical (pendiente infinita, es decir x_diff == 0)
				for(size_t y = 0; y < _img.rows(); y++){
					//x0 == x1
					std::pair<size_t, size_t> celda = cell_from_pixel(y, x0);
					size_t posicion_celda = celda.first*_n + celda.second;
					D(i, posicion_celda) += 1;
					velocidades_rayo_actual[posicion_celda] += _img(y, x0);
				}
			}
			else{

				for(size_t x = 0; x < _img.cols(); x++){ // _img es cuadrada, es lo mismo rows y cols
					
					int y = round( ((x-x0)/x_diff)*y_diff + y0 );
					if (y < 0  || y >= (long int)_img.rows()) continue; //No sirve ningún pixel por fuera de la discretización.
					
					std::pair<size_t, size_t> celda = cell_from_pixel(y,x);
					if(celda.first > _n || celda.second > _n) continue; //No sirve ninguna celda por fuera de la discretización
					
					size_t posicion_celda = celda.first*_n + celda.second;
					D(i, posicion_celda) += 1;
					velocidades_rayo_actual[posicion_celda] += _img(y,x);
				}
			}

			//Calcular t y reiniciar vector de velocidades para el siguiente rayo
			for(size_t k = 0; k < velocidades_rayo_actual.size(); k++){
				t(i,0) += (velocidades_rayo_actual[k]);
				velocidades_rayo_actual[k] = 0;
			}

			
		}
		return resultado;
	}


	/*
	Dado un pixel, devuelve la celda a la que pertenece en la discretización. Es 0-indexed.
	pre: i y j van de 0 a _img.rows() - 1. _img es una matriz cuadrada.
	*/
	std::pair<size_t, size_t> cell_from_pixel(size_t i, size_t j){
		unsigned ppc = _img.rows() / _n; //Pixels Por Celda
		size_t fila	   = (size_t)(i/ppc); 
		size_t columna = (size_t)(j/ppc);
		return std::make_pair(fila, columna);
	}

	/*
	Carga una imagen ppm desde un archivo en la discretizacion
	*/
	Matrix load_matrix(const std::string& filename ){
		
		uchar* data = nullptr;
		int width = 0, height = 0;
		PPM_LOADER_PIXEL_TYPE pt = PPM_LOADER_PIXEL_TYPE_INVALID;
		bool ret = LoadPPMFile(&data, &width, &height, &pt, filename.c_str());
		if(!ret || width == 0|| height == 0) throw std::runtime_error("no se pudo abrir la imagen");

		if( pt != PPM_LOADER_PIXEL_TYPE_GRAY_8B && pt != PPM_LOADER_PIXEL_TYPE_GRAY_16B){
			throw std::runtime_error("Profundida de bits no soportada");
		}

		unsigned short* s_data = (unsigned short*) data;

		Matrix loaded_matrix(height, width);
		for (int h = 0; h < height; ++h){
			for (int w = 0; w < width; ++w){
				int colIndex = h * width + w;
				switch(pt){
				
					case PPM_LOADER_PIXEL_TYPE_GRAY_8B:
						loaded_matrix(h, w) = data[colIndex];		
					break;

					case PPM_LOADER_PIXEL_TYPE_GRAY_16B:
						loaded_matrix(h, w) = s_data[colIndex];		
					break;
				}

			}
		}
		return loaded_matrix;
	}

};

#endif