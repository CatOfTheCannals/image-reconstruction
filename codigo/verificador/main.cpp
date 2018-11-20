
#include "matrix.h"
#include "matrix_algos.h"
#include "../ppmloader/ppmloader.h"
#include "discretizador.h"
#include "geometrias.h"

#include <iostream>
#include <chrono>
using namespace std;

std::vector<std::string> str_tokenize(const std::string& str, const std::string& characters){

	std::vector<std::string> tokens;
	if( str.size() != 0){
		size_t initial_pos = 0;
		size_t ocurrence = str.find_first_of(characters, initial_pos);
		if(ocurrence == std::string::npos) tokens.push_back(str);
		while(ocurrence != std::string::npos){
			std::string s = str.substr(initial_pos, ocurrence - initial_pos);
			if (s.size() > 0) tokens.push_back(s);
			initial_pos = ocurrence + 1;
			ocurrence = str.find_first_of(characters, initial_pos);
		}
		if(initial_pos < str.size()) tokens.push_back(str.substr(initial_pos));
	}

	return tokens;
}

void print_help(char* filename){
	cout << endl;
	cout << "uso: " << endl;
	cout << filename << " -i <input> -o <output> [-r <ruido> -gX [params]]" << endl;
	cout << endl;
	cout << "donde: "		<< endl;
	cout << "<input>\t"	<< "Es la imagen de entrada. Debe ser cuadrada." << endl;
	cout << "<output>\t"	<< "Es la imagen de salida. Es cuadrada con una resolucion de <n>x<n>." << endl;
	cout << "<ruido>\t" << "nivel de ruido. Agrega ruido aleatorio entre 0 y <ruido>. Si no especifica se usa 0.0." << endl;
	cout << "gX determina la geometria a usar. [params] son los parametros de la geometria." << endl;
	cout << endl;
}


int main(int argc, char** argv){
	
	if(argc == 1 || argc % 2 != 1){
		cout << "Cantidad de parámetros inválida" << endl;
		print_help(argv[0]);
		return 1;
	}

	string input_img;
	string output_img;
	bool debug = false;
	vector<pair<string, string>> geometrias;
	double r(0); //nivel de ruido

	for(int i = 1; (i+1) < argc ; i+=2 ){
		string param(argv[i]); 
		
		if(param == "-i"){
			input_img = argv[i+1];
		}
		else if(param == "-o"){
			output_img = argv[i+1];
		}
		else if(param == "-r"){
			r = stod(argv[i+1]);	
		}
		else if(param == "-d"){
			debug = true;
		}
		else if(param.substr(0, 2) == "-g"){
			geometrias.push_back(make_pair(param, argv[i+1]));
		}
		else{
			cout << r << ": parametro desconocido\n";
			exit(2);
		}
	}

	if(input_img.size() == 0 || output_img.size() == 0){
		cout << "parámetros inválidos\n";
		print_help(argv[0]);
		exit(1);
	}


	try{

		Discretizador discretizacion(input_img, 1, r);

		std::chrono::system_clock::time_point begin = std::chrono::system_clock::now();
		std::vector<std::string> g5_params = {"", ""};
		std::vector<std::string> g14_params = {"", ""};
		std::vector<std::string> g15_params = {"", ""};
		for( auto geo : geometrias ){

			string num_geo = geo.first;
			if(num_geo == "-g1" ){
				geo_horizontales h(stoi(geo.second));
				h.aplicar(discretizacion);
			} 
			else if(num_geo == "-g2" ){
				geo_verticales v(stoi(geo.second));
				v.aplicar(discretizacion);
			}
			else if(num_geo == "-g3"){
				geo_diagonales d(stoi(geo.second));
				d.aplicar(discretizacion);
			}
			else if(num_geo == "-g4"){
				geo_random r(stoi(geo.second));
				r.aplicar(discretizacion);
			}
			else if(num_geo == "-g10"){
				geo_horizontalesPorCantidad hc(stoi(geo.second));
				hc.aplicar(discretizacion);
			}
			else if(num_geo == "-g11"){
				geo_verticalesPorCantidad vc(stoi(geo.second));
				vc.aplicar(discretizacion);
			}
			else if(num_geo == "-g12"){
				geo_diagonalesPorCantidad dc(stoi(geo.second));
				dc.aplicar(discretizacion);
			}
			else if(num_geo == "-g13"){
				geo_centralPorCantidad cc(stoi(geo.second));
				cc.aplicar(discretizacion);
			}
			else if(num_geo == "-g5"){
				auto params = str_tokenize(geo.second, ",");
				geo_punto_fijo_random pf(stoi(params[0]), stoi(params[1])); //ejemplo ./tp3 [...] -g5 5,6 -> llama pf(5, 6)
				pf.aplicar(discretizacion);
			}
			else if(num_geo == "-g6"){
				geo_central c(stoi(geo.second));
				c.aplicar(discretizacion);
			}
			else if(num_geo == "-g7"){
				auto params = str_tokenize(geo.second, ",");
				geo_focosRegularesIzquierdaDerecha fr(stoi(params[0]), stoi(params[1]));
				fr.aplicar(discretizacion);
			}
			else if(num_geo == "-g8"){
				auto params = str_tokenize(geo.second, ",");
				geo_focosRegularesArribaAbajo fr(stoi(params[0]), stoi(params[1]));
				fr.aplicar(discretizacion);
			}
			else if (num_geo == "-g9"){
				auto params = str_tokenize(geo.second, ",");
				discretizacion.trazar_rayo(stoi(params[1]), stoi(params[0]), stoi(params[3]), stoi(params[2]));
			}
			else if(num_geo == "-g10"){
                geo_horizontalesPorCantidad hc(stoi(geo.second));
                hc.aplicar(discretizacion);
            }
            else if(num_geo == "-g11"){
                geo_verticalesPorCantidad vc(stoi(geo.second));
                vc.aplicar(discretizacion);
            }
                else if(num_geo == "-g12"){
                geo_diagonalesPorCantidad dc(stoi(geo.second));
                std::cout << "debug: rays applied" << std::endl;
                dc.aplicar(discretizacion);
                std::cout << "debug: geo applied" << std::endl;
            }
            else if(num_geo == "-g13"){
                geo_centralPorCantidad cc(stoi(geo.second));
                cc.aplicar(discretizacion);
            }
            else if(num_geo == "-g14"){
            	auto params = str_tokenize(geo.second, ",");
            	geo_focosRegularesIzquierdaDerechaFijos  fridf(stoi(params[0]), stoi(params[1]));
            	fridf.aplicar(discretizacion);
            }
            else if(num_geo == "-g14a"){
				g14_params[0] = geo.second;
				//auto params = str_tokenize(geo.second, ",");
				if(g14_params[1] != ""){
					geo_focosRegularesIzquierdaDerechaFijos pf(stoi(g14_params[0]), stoi(g14_params[1])); //ejemplo ./tp3 [...] -g5 5,6 -> llama pf(5, 6)
					pf.aplicar(discretizacion);
					g14_params[0] = "";
					g14_params[1] = "";
				}
			}
			else if(num_geo == "-g14b"){
				g14_params[1] = geo.second;
				//auto params = str_tokenize(geo.second, ",");
				if(g14_params[0] != ""){
					geo_focosRegularesIzquierdaDerechaFijos pf(stoi(g14_params[0]), stoi(g14_params[1])); //ejemplo ./tp3 [...] -g5 5,6 -> llama pf(5, 6)
					pf.aplicar(discretizacion);
					g14_params[0] = "";
					g14_params[1] = "";
				}	
			}
            else if(num_geo == "-g15"){
            	auto params = str_tokenize(geo.second, ",");
            	geo_focosRegularesArribaAbajoFijos  fridf(stoi(params[0]), stoi(params[1]));
            	fridf.aplicar(discretizacion);
            }
            else if(num_geo == "-g15a"){
				g15_params[0] = geo.second;
				//auto params = str_tokenize(geo.second, ",");
				if(g15_params[1] != ""){
					geo_focosRegularesArribaAbajoFijos pf(stoi(g15_params[0]), stoi(g15_params[1])); //ejemplo ./tp3 [...] -g5 5,6 -> llama pf(5, 6)
					pf.aplicar(discretizacion);
					g15_params[0] = "";
					g15_params[1] = "";
				}
			}
			else if(num_geo == "-g15b"){
				g15_params[1] = geo.second;
				//auto params = str_tokenize(geo.second, ",");
				if(g15_params[0] != ""){
					geo_focosRegularesArribaAbajoFijos pf(stoi(g15_params[0]), stoi(g15_params[1])); //ejemplo ./tp3 [...] -g5 5,6 -> llama pf(5, 6)
					pf.aplicar(discretizacion);
					g15_params[0] = "";
					g15_params[1] = "";
				}	
			}
		}

	  	std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
	  	double tiempo_rayos = std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

		discretizacion.verificar_trazado(output_img, debug);

		
	}
	catch(runtime_error e){
		cout << e.what() << endl;
	}

	return 0;
}

