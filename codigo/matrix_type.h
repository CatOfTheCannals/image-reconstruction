#ifndef MATRIX_TYPE_H
#define MATRIX_TYPE_H

#include <assert.h>
#include <iostream>
#include <memory>
#include <cstring>
#include <vector>
#include <limits>
#include <math.h>
#include <fstream>

#ifndef BASE_INDEX
	#define BASE_INDEX 0
#endif

template <typename T>
class matrix_type{

private:

	size_t _rows;
	size_t _cols;
	std::string _name;
	std::unique_ptr<T[]> _data;


public:

	typedef T value_type;
	const T epsilon;

	matrix_type(): matrix_type(1,1){
		insert(BASE_INDEX, BASE_INDEX, value_type(1));
	}

	matrix_type(size_t i, size_t j):
		_rows(i),
		_cols(j),
		_name("Matrix"),
		_data(new T[i*j]),
		epsilon(0.000004){

		assert(i != 0 && j!= 0);
	}

	matrix_type(const matrix_type<T>& other):
		matrix_type(other._rows, other._cols){

		_name = other._name;
		for(size_t i = 0; i < _rows*_cols; i++)
			_data[i] = other._data[i];
	}

	matrix_type(const std::string& file_name):
		_name("Matrix"),
		epsilon(0.000004){
		if(!load(file_name)){
			_rows = 1;
			_cols = 1;
			_data.reset(new T[1]);
			throw false;
		}
	}

	matrix_type<T>& operator=(matrix_type other){
		swap(*this, other);
		return *this;
	}

	inline
	size_t rows() const{
		return _rows;
	}

	inline
	size_t cols() const{
		return _cols;
	}

	void set_name(const std::string& name){
		_name = name;
	}
	
	const std::string& name() const{ 
		return _name; 
	}

	inline 
	void insert(size_t i, size_t j, T element){

		assert(i - BASE_INDEX < _rows && j - BASE_INDEX < _cols);
		_data[_cols*(i-BASE_INDEX) + (j - BASE_INDEX)] = element;
	}

	inline 
	T& operator()(size_t i, size_t j){

		assert(i - BASE_INDEX < _rows && j - BASE_INDEX < _cols);
		return _data[_cols*(i-BASE_INDEX) + (j - BASE_INDEX)];
	}

	inline 
	const T& operator()(size_t i, size_t j) const{

		assert(i - BASE_INDEX < _rows && j - BASE_INDEX < _cols);
		return _data[_cols*(i-BASE_INDEX) + (j - BASE_INDEX)];
	}


	bool store(const std::string& file_name){

		try{
			std::ofstream output(file_name);
			output << rows() << std::endl;
			output << cols() << std::endl;
			for(size_t i = 0; i < rows()*cols(); i++) output << _data[i] << std::endl;
			output.close();
		}
		catch(std::exception e){
			std::cout << "ERROR: " << e.what() << std::endl;
			return false;
		}
		return true;
	}

	bool load(const std::string file_name){
		try{
			std::ifstream input(file_name);
			if (!input.good()) return false;
			input >> _rows;
			input >> _cols;
			_data.reset(new T[_rows*_cols]);
			for(size_t i = 0; i < rows()*cols(); i++) input >> _data[i];
			input.close();
		}
		catch(std::exception e){
			std::cout << "ERROR: " << e.what() << std::endl;
			return false;	
		}
		return true;
	}

	matrix_type<T> mt_times_m() const{
		//X.transpose()*X
		matrix_type result = matrix_type(this->cols(), this->cols());

		int index = 0;
		for (int i = 0; i < this->cols(); i++) {
			double* i_column = new double[this->rows()];
			for (int k = 0; k < this->rows(); k++) {
				i_column[k] = (*this)(k, i);
			}

			for (int j = 0; j < this->cols(); j++) {

				double temp = 0;
				//actual multiplication
				for (int k = 0; k < this->rows(); k++) {
					temp += i_column[k] * (*this)(k, j);
				}

				result._data[index] = temp;
				index++;
			}
			delete[] i_column;

		}
		return result;
	}

	friend
	void swap(matrix_type& first, matrix_type& second){
		using std::swap;
		swap(first._data, second._data);
		swap(first._rows, second._rows);
		swap(first._cols, second._cols);
		swap(first._name, second._name);
	}

	friend
	void eliminacion_gausseana(matrix_type& A){

		size_t n = A.rows() - 1 + BASE_INDEX;
		size_t m = A.cols() - 1 + BASE_INDEX;

		for(size_t j = BASE_INDEX; j <= n - 1; j++){
			for(size_t i = j + 1; i <= n; i++){
				if ( fabs(A(i,j)) < std::numeric_limits<T>::min()*10) continue;
				auto M = A(i,j)/A(j,j);
				for(size_t k = j; k <= m; k++){
					A.insert(i, k , A(i,k) - M*A(j,k));
				}
			}
		}
	}



};

#endif
