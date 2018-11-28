#ifndef MATRIX_ALGOS_CPP
#define MATRIX_ALGOS_CPP

#include "matrix_algos.h"

/* Imprime una matriz */
std::ostream& operator<<(std::ostream& os, const Matrix& M){

	for(int i = 0; i < M.rows() ; i++){
        for(int j = 0; j < M.cols(); j++){
            os << M(i,j) << " ";
        }
        os << std::endl;
    }

	return os;
}

/*Multiplicación de Matrix por escalar */
Matrix scalar_mult(const typename Matrix::value_type& scalar, const Matrix& A){

	Matrix C(A.rows(), A.cols());
	for(size_t i = BASE_INDEX; i < A.rows() + BASE_INDEX; i++)
		for(size_t j = BASE_INDEX; j < A.cols() + BASE_INDEX; j++)
			C.insert(i,j ,scalar*A(i,j));

	return C;
}

/*Suma de matrices */
Matrix operator+(const Matrix& A, const Matrix& B){	

	assert(A.rows() == B.rows() && A.cols() == B.cols());

	Matrix C(A.rows(), A.cols());
	for(size_t i = BASE_INDEX; i < A.rows() + BASE_INDEX; i++)
		for(size_t j = BASE_INDEX; j < A.cols() + BASE_INDEX; j++)
			C.insert(i,j, A(i,j) + B(i,j));

	return C;

}

Matrix operator-(const Matrix& A, const Matrix& B){	
	
	assert(A.rows() == B.rows() && A.cols() == B.cols());

	Matrix C(A.rows(), A.cols());
	for(size_t i = BASE_INDEX; i < A.rows() + BASE_INDEX; i++)
		for(size_t j = BASE_INDEX; j < A.cols() + BASE_INDEX; j++)
			C.insert(i, j, A(i,j) - B(i,j));

	return C;

}

bool operator==(const Matrix& A, const Matrix& B){

	assert(A.rows() == B.rows() && A.cols() == B.cols());


	for(size_t i = BASE_INDEX; i < A.rows() + BASE_INDEX; i++)
		for(size_t j = BASE_INDEX; j < A.cols() + BASE_INDEX; j++){
			if( fabs(A(i,j) - B(i,j)) >  A.epsilon){ //todas las matrices tienen el mismo epsilon
				return false;
			}
		}

	return true;

}

bool operator!=(const Matrix& A, const Matrix& B){

	return !(A == B);
}

/*Multiplicacion de 2 submatrices principales de A y B */
//l: filas de submatriz de A
//m: columnas de submatriz de A y filas de submatriz de B
//n: columnas de submatriz de B
Matrix partial_matrix_mult(const Matrix& A, const Matrix& B, size_t l, size_t m, size_t n){
	
	l = std::min(l, A.rows());
	m = std::min(m, A.cols());
	n = std::min(n, B.cols());

	assert(m <= B.rows());

	Matrix C = get_zero_initialized(l, n);

	for(size_t i = BASE_INDEX; i < l + BASE_INDEX; i++){
		for(size_t j = BASE_INDEX; j < m + BASE_INDEX; j++){
			for(size_t k = BASE_INDEX; k < n + BASE_INDEX; k++){
				C.insert(i, k, C(i,k) + A(i, j)*B(j,k));
			}

		}
	}

	return C;
}

Matrix operator*(const Matrix& A, const Matrix& B){
	
	assert(A.cols() == B.rows());
	Matrix C(A.rows(), B.cols());

	for(size_t i = BASE_INDEX; i < A.rows() + BASE_INDEX; i++){
		for(size_t j = BASE_INDEX; j < B.cols() + BASE_INDEX; j++){
            double acum = 0;
            for(int k = BASE_INDEX; k < A.cols() + BASE_INDEX ; k++){
			    acum+= A(i,k)*B(k,j);
            }
            C(i,j) = acum;
		}
	}

	return C;
}

Matrix operator/(const Matrix& A, const double& scalar){
    Matrix C(A.rows(), A.cols());
    for(size_t i = BASE_INDEX; i < A.rows() + BASE_INDEX; i++)
        for(size_t j = BASE_INDEX; j < A.cols() + BASE_INDEX; j++)
            C.insert(i,j ,A(i,j)/scalar);

    return C;
}

//B es una matrix cuadrada diagonal
Matrix diagonal_mult(const Matrix& A, const Matrix& B){

	assert(A.rows() == B.rows() && A.cols() == B.cols());

	
	Matrix C(A);
	size_t n = C.rows() - 1 + BASE_INDEX;
	size_t m = C.cols() - 1 + BASE_INDEX;

	std::vector<Matrix::value_type> d(B.rows());
	for(size_t i = 0; i < d.size(); i++) 
		d[i] = B(i + BASE_INDEX, i + BASE_INDEX);

	for(size_t i = BASE_INDEX; i <= n; i++){
		for(size_t j = BASE_INDEX; j <= m; j++){
			C.insert(i, j, C(i,j)*d[j - BASE_INDEX]);
		}
	}

	return C;
}


Matrix get_initialized(size_t n, size_t m, Matrix::value_type element){
	Matrix C(n, m);
	for(size_t i = BASE_INDEX; i < n + BASE_INDEX; i++)
		for(size_t j = BASE_INDEX; j < m + BASE_INDEX; j++)
			C.insert(i, j ,element);
	return C;
}

Matrix get_zero_initialized(size_t n, size_t m){
	return get_initialized(n, m, Matrix::value_type(0));
}

Matrix get_zero_initialized(size_t n){
	return get_zero_initialized(n, n);
}

/*Devuelve una matriz identidad */
Matrix get_identity(size_t n){
	Matrix::value_type one = 1;
	Matrix C = get_zero_initialized(n);

	for(size_t i = BASE_INDEX; i < n + BASE_INDEX; i++)
		C.insert(i,i ,one);
	return C;
}


//Requiere que B sea un vector columna con la misma cantidad de filas que A
Matrix augmentar(const Matrix& A, const Matrix& B){

	assert(B.cols() == 1 && A.rows() == B.rows());
	Matrix C(A.rows(), A.cols() + 1);

	for(size_t i = BASE_INDEX; i < A.rows() + BASE_INDEX; i++){
		for(size_t j = BASE_INDEX; j < A.cols() + BASE_INDEX; j++){
			C.insert(i, j, A(i,j));
		}
	}

	for(size_t i = BASE_INDEX; i < B.rows() + BASE_INDEX; i++){
		C.insert(i, A.cols() + BASE_INDEX, B(i, BASE_INDEX));
	}
	return C;
}

//Requiere que A sea la matrix resultante de aplicar eliminación gausseana sobre una matrix augmentada
Matrix backward_substitution(const Matrix& A){

	size_t n = A.rows() - 1 + BASE_INDEX;
	size_t m = A.cols() - 1 + BASE_INDEX;

	Matrix X(A.rows(), 1);

	X.insert(n, BASE_INDEX, A(n, m)/A(n, m-1));

	for(size_t i_aux = BASE_INDEX; i_aux <= n - 1; i_aux++){
		size_t i = n - i_aux;

		//Calcular x_i
		//x_i = b_i - sum(a_ij con j > i) / a_ii
		Matrix::value_type sum(0), c(0);
		for(size_t j = i + 1; j <= m - 1; j++){
			Matrix::value_type y = A(i, j)*X(j, BASE_INDEX) - c;
			Matrix::value_type t = sum + y;
			c = (t - sum) - y;
			sum = t;
		}

		X.insert(i, BASE_INDEX, (A(i, m) - sum)/A(i,i));
	}
	return X;
}

Matrix resolver_sistema(const Matrix& A, const Matrix& b){

	assert(A.cols() == b.rows() && b.cols() == 1);
	Matrix tmp = augmentar(A, b);
	eliminacion_gausseana(tmp);
	return backward_substitution(tmp);
}

std::tuple<Matrix, Matrix, double, int> generar_svd(const Matrix& A){
    Matrix AtA = A.mt_times_m();
    auto v_s_k = calcular_autovectores(AtA);
    Matrix v = std::get<0>(v_s_k);
    Matrix s = std::get<1>(v_s_k);
    double k = std::get<2>(v_s_k);

	matrix_stats(v, "v");
	matrix_stats(s, "s");

	Matrix sq = sqrt_to_all_elems(s);
	double numero_de_condicion = sq(BASE_INDEX, BASE_INDEX) / sq(BASE_INDEX + sq.rows()-1, BASE_INDEX);

    // v tiene k columnas
    Matrix u_s = A * v;

    Matrix ut = apply_inverse_sigma(sq, trasponer(u_s), k);
	std::cout << "ut built" << std::endl;

    Matrix sInv_ut = apply_inverse_sigma(sq, ut, k);
	std::cout << "sInv_ut built" << std::endl;

    return std::make_tuple(v, sInv_ut, numero_de_condicion, k);
}

Matrix resolver_sistema_con_svd(const Matrix& v, const Matrix& sInv_ut, const Matrix& b){
	Matrix v_sInv_ut_b = v * (sInv_ut * b);
	std::cout << "v_sInv_ut_b built" << std::endl;
    return v_sInv_ut_b;
}

Matrix sqrt_to_all_elems(const Matrix& A) {
	Matrix res(A);
	for (int  i = 0;  i < res.rows(); ++ i) {
		for (int j = 0; j < res.cols(); ++j) {
			res.insert(i,j, sqrt(res(i,j)) );
		}
	}
	return res;
}

void matrix_stats(const Matrix& A, const std::string name) {
    double epsilon = 0.0001;
    int nans = 0;
    int zeros= 0;
    int negatives = 0;

    std::cout << "debug: __matrix_stats__" << name << std::endl;

    std::cout << "debug: rows " << A.rows() << std::endl;
    std::cout << "debug: cols " << A.cols() << std::endl;

    for (int  i = 0;  i < A.rows(); ++ i) {
        for (int j = 0; j < A.cols(); ++j) {
            if(isnan(A(i,j))) {
                nans ++;
            } else if(fabs(A(i,j)) < epsilon ) {
                zeros ++;
            } else if(A(i,j) < 0) {
                negatives ++;
            } else {
                // std::cout << "debug: s " << i << " = " << A(i,0) << std::endl;
            }
        }
    }
    std::cout << "debug: nans " << nans << std::endl;
    std::cout << "debug: zeros " << zeros << std::endl;
    std::cout << "debug: negatives " << negatives << std::endl;

}

bool has_nan(const Matrix& A) {
    bool res = false;
    for (int  i = 0;  i < A.rows(); ++ i) {
        for (int j = 0; j < A.cols(); ++j) {
            if(isnan(A(i,j))) {
                res = true;
                break;
            }
        }
    }
    return res;
}

bool is_zero(const Matrix& A) {
    bool res = true;
    double epsilon = 0.0001;

    for (int  i = 0;  i < A.rows(); ++ i) {
        for (int j = 0; j < A.cols(); ++j) {
            if(fabs(A(i,j)) > epsilon ) {
                res = false;
                break;
            }
        }
    }
    return res;
}


Matrix apply_inverse_sigma(const Matrix& s, const Matrix& A, int rows) {

    Matrix res = subMatrix(A, 0, rows, 0, A.cols()-1);
	std::cout << "submatrix built " << std::endl;

	matrix_stats(res,"res");
	matrix_stats(s,"s");

    // for i, singular_value in enumerate(s):
    for (int  i = 0;  i < s.rows(); ++ i) {
        for (int j = 0; j < res.cols(); ++j) {

			std::cout << std::endl;
			std::cout << "res.rows()" << res.rows() << "r es.cols()" << res.cols() << std::endl;
			std::cout << "s.rows()" << s.rows() << " s.cols()" << s.cols() << std::endl;
			std::cout << "i" << i << " j" << j << std::endl;

			double num = res(i,j);
			std::cout << "num calculated " << std::endl;

			double den = s(i, 0);
			std::cout << "den calculated " << std::endl;

			std::cout << "den" << den << std::endl;

			double val = num / den;
			std::cout << "val calculated " << std::endl;

			if(isnan(val)) {
				std::cout << "debug: nan appeared at " << i << ", " << j << std::endl;
				std::cout << "debug: previous val " << res(i,j) << std::endl;
				std::cout << "debug: sigma " << s(i,0) << std::endl;
				std::cout << std::endl;
			}
			res.insert(i,j, val);
			std::cout << "val inserted " << std::endl;
        }
    }
    return res;
}

Matrix subMatrix(const Matrix& A, int i1, int i2, int j1, int j2){
    /*
     * returns the matrix between (i1, i2) rows and (j1, j2) cols
     * indexes i2 and j2 will be part of the answer
     */

    assert(i1 <= i2 && j1 <= j2);
    assert(-1 < i1 && -1 < j2);
    assert(i2 < A.rows() && j2 < A.cols());

    int res_rows = i2 - i1 + 1;
    int res_cols = j2 - j1 + 1;

    Matrix res(res_rows, res_cols);
    for(int i = BASE_INDEX; i < res.rows() + BASE_INDEX; i++){
        for(int j = BASE_INDEX; j < res.cols() + BASE_INDEX; j++){
            res.insert(i,j, A(i + i1, j1 + j));
        }
    }
    return res;
}

bool is_relevant(double d){
    double epsilon = 1e-5;
    if(d < 0){ return d < epsilon*(-1); }
    else{ return d > epsilon; }
}

Matrix producto_externo(double lambda, Matrix& v){
    assert(v.cols()==1);
    Matrix res = Matrix(v.rows(),v.rows());
    for(int i = 0 ; i < v.rows() ; i++){
        for(int j = 0 ; j < v.rows() ; j++){
            res.insert(i, j, lambda*v(i,BASE_INDEX)*v(j,BASE_INDEX));
        }
    }
    return res;
}


//Calcula los primeros k autovectores de B
std::tuple<Matrix, Matrix, double> calcular_autovectores(Matrix B){
	assert(B.rows() == B.cols());
	size_t n = B.rows();

	//El vector inicial se calcula random.
	//Ajustamos iteraciones según nuestro criterio.
	Matrix v(n, 1);
	Matrix autovalores(n, 1);
	Matrix v_anterior(n, 1);
	Matrix C(n, n);
	int iteraciones = 100;
	int MAX_RETRIES = 5;
    int retries = 0;
	int k = n-1;
	//calcular n autovectores
	for(size_t autovector_actual = BASE_INDEX; autovector_actual < n + BASE_INDEX; autovector_actual++){
		for(size_t i = BASE_INDEX; i < n + BASE_INDEX; i++){
			// v.insert(i, BASE_INDEX, (1 - 2*(i%2)) );
			v.insert(i, BASE_INDEX, rand()%100 );
		}

		for(int i = 0; i < iteraciones ; i++){
			v_anterior = v;
			v = B*v;

			double norma = norma_2(v);
			double inv_norm = 1/norma;
			v = scalar_mult(inv_norm, v);
		}

		Matrix v_t = trasponer(v);

		for(size_t i = BASE_INDEX; i < n + BASE_INDEX; i++){
			C.insert(i, autovector_actual, v(i, BASE_INDEX));
		}

		Matrix numer = B*v;
		double autovalor_asociado = 0;
		if(retries==0){
            autovalor_asociado = (v_t*numer)(BASE_INDEX,BASE_INDEX);
            double denom = (v_t*v)(BASE_INDEX,BASE_INDEX); //hopefully, about 1.
            autovalor_asociado = autovalor_asociado/denom;
		} else {
			double norm_coef = 10;
			for(int retry = 0; retry < retries; retry++){norm_coef*=10;}
			Matrix v_t_primo = scalar_mult(norm_coef, v_t);
		    Matrix numer_primo = scalar_mult(norm_coef, numer);
            autovalor_asociado = (v_t_primo*numer_primo)(BASE_INDEX,BASE_INDEX);
            autovalor_asociado /= (norm_coef*norm_coef);
		}

		double autovalor_segundo = numer(BASE_INDEX+1, BASE_INDEX)/v(BASE_INDEX+1, BASE_INDEX);

		/*
		    get_eigenvalue(numer, v);
		*/
		autovalores.insert(autovector_actual, BASE_INDEX, autovalor_asociado);

		std::cout << v(BASE_INDEX,BASE_INDEX) << ", " << v(BASE_INDEX+1, BASE_INDEX) << std::endl ;
		std::cout << "autovalores : " << autovalor_asociado <<", " << autovalor_segundo << std::endl ;
		std::cout << numer(BASE_INDEX,BASE_INDEX)<< ", " << numer(BASE_INDEX+1, BASE_INDEX) << std::endl ;

		bool ruined = false;
		if(autovalor_asociado < 0 || autovalor_segundo < 0){
            ruined = true;
			std::cout << autovalor_asociado << ", " << autovalor_segundo << std::endl ;
			std::cout << "We fucked up. Trying again." << std::endl ;
		}
		if(!ruined){
			Matrix extern_prod = producto_externo(autovalor_asociado, v);
		    B = B - extern_prod;
		    retries = 0;
		} else {
			retries++;
			std::cout << "Current number of retries: " << retries << std::endl ;
			if(retries > MAX_RETRIES){
				std::cout << "We fucked up too much. Suicide is the only option now" << std::endl;
				k = autovector_actual - 1;
				break;

			} else {
				autovector_actual--;
			}
			
		}
		
	}

    // achicar C
	Matrix useful_autovecs = subMatrix(C, 0, n-1, 0, k-1);
    // achicar autovalores
	Matrix useful_autovals = subMatrix(autovalores, 0, k, 0, 0);

	return std::make_tuple(useful_autovecs, useful_autovals, k);
}

double get_eigenvalue(const Matrix& numer, const Matrix& v){
	double autovalor_asociado = 0;
	bool escape = false;
	for(int val = 0; val < v.rows() && escape==false; val++){
		double denom = v(val,0);
		if( is_relevant(denom) ){
			autovalor_asociado = numer(val,0)/denom;
		    escape = true;
		}
	}
	return autovalor_asociado;
}

Matrix random(int rows, int cols) {
    Matrix res(rows, cols);
    srand (time(NULL));
    for (std::size_t i = 0; i < rows; i++) {
        for (std::size_t j = 0; j < cols; j++) {
            res.insert(i, j, rand() % 100);
        }
    }
    return res;
}

Matrix trasponer(const Matrix& A){

	Matrix C(A.cols(), A.rows());

	for(size_t i = BASE_INDEX; i < A.rows() + BASE_INDEX; i++)
		for(size_t j = BASE_INDEX; j < A.cols() + BASE_INDEX; j++)
			C.insert(j, i, A(i,j));

	return C;
}

Matrix::value_type norma_2(const Matrix& v){
	assert(v.cols() == 1);

	size_t n = v.rows();	
	double res = 0;

	for(size_t i = BASE_INDEX; i < n + BASE_INDEX; i++){
		double y = v(i, BASE_INDEX);
		res += y*y;
	}
	return sqrt(res);
}

inline
Matrix::value_type norma_inf(const Matrix& v){
	assert(v.cols() == 1);

	typedef Matrix::value_type decimal;

	size_t n = v.rows();	
	decimal norma_inf(0);
	decimal c(0);
	for(size_t i = BASE_INDEX; i < n + BASE_INDEX; i++){
		decimal y = std::abs(v(i, BASE_INDEX)) - c;
		decimal t = norma_inf + y;
		c = (t - norma_inf) - y;
		norma_inf = t;
	}

	return norma_inf;

}

#endif
