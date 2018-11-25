#include "../codigo/matrix.h"
#include "../codigo/matrix_algos.h"
//#include "../codigo/discretizador.h"
//#include "../codigo/geometrias.h"
//#include "../codigo/ppmloader/ppmloader.h"

#include <iostream>
#include <chrono>
using namespace std;

void test_norma_2(){
	cout << "test norma 2 ..." << endl;
	Matrix v = Matrix(5,1);
	for(int i = 0; i< 5; i++){v.insert(i, BASE_INDEX, i+1);}
	cout << v << endl ;
    double size = norma_2(v);
    cout <<"norma 2: "<< size << " ~ 7.41..." << endl ;
    cout << "fin test norma 2 ..." << endl;

}

void test_prod_ext(){
	cout << "test producto externo ..." << endl;
	Matrix v = Matrix(5,1);
	for(int i = 0; i< 5; i++){v.insert(i, BASE_INDEX, i+1);}
	cout << v << endl ;
    double d = 0.0023;
    Matrix ext = producto_externo(d, v);
    cout <<"v*v_t : \n"<< ext << endl ;
    cout << "fin test producto externo ..." << endl;

}

void test_prod_normal(){
	cout << "test producto de matrices ..." << endl;
	Matrix m = Matrix(5,5);
	for(int i = 0; i< 5; i++){
		for(int j = 0; j< 5; j++){
			m.insert(i,j, i+j);
		}
	}
	cout << m << endl ;
    Matrix m_2 = m*m;
    cout <<"m * m : \n"<< m_2 << endl ;
    cout << "fin test producto de matrices ..." << endl;

}

void test_autovalores(){
	cout << "test autovalores ..." << endl;
	int DIM = 10;
    Matrix testy = Matrix(DIM,DIM);

    for(int i = 0; i< DIM; i++){testy.insert(i, i, DIM+1-i);}
    cout << testy << endl;
	tuple<Matrix, Matrix> ev_ev = calcular_autovectores(testy, testy.rows());

    Matrix autovalores = get<1>(ev_ev);
    cout << "autovalores rows: " << autovalores.rows() << endl;
    cout << "autovalores: \n " << autovalores << endl;
    Matrix autovectores = get<0>(ev_ev);
    //cout << "autovectores: "<< autovectores << endl;
    cout << "fin test autovalores ..." << endl;
}

int main(){
	test_norma_2(); //TEST OK 
	test_prod_normal();
	test_prod_ext();
    test_autovalores();

	return 0;
}