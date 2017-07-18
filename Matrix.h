#ifndef Matrix_h_
#define Matrix_h_

#include <cstdlib> // For exit()
#include <iostream> // For cerr
#include <vector>

using std::vector;

template<typename DataType>
class Matrix {
	typedef vector<DataType> RowType;

public:
	size_t number_of_rows;
	size_t number_of_columns;

	Matrix();
	Matrix(size_t number_of_rows, size_t inColumns, const DataType& value =
			DataType());
	Matrix(const Matrix& matrix);

	~Matrix(void) {
	}

	vector<size_t> size() {
		vector<size_t> size(2);
		size.at(0) = number_of_rows;
		size.at(1) = number_of_columns;
		return size;
	}

	bool empty() const {
		return number_of_rows == 0;
	}

	void fill(const DataType& value);
	void clear();

	void resize(size_t inRows, size_t inColumns, const DataType& inVal =
			DataType());

	const RowType& operator[](size_t index) const;
	RowType& operator[](size_t index);

	DataType& at(size_t row_index, size_t column_index);
	DataType& at(size_t row_index, size_t column_index) const;

	const Matrix& operator+(const Matrix& matrix);
	const Matrix& operator+=(const Matrix& matrix);

	void SwapRows(size_t row_1, size_t row_2);
	void SwapColumns(size_t column_1, size_t column_2);

	void swap(Matrix& matrix);
	void symmetrize();

private:
	vector<RowType> rows;
};

/**
 * STP: This default constructor is not necessary because all of the members
 * will be initialized with their default constructors anyway.
 *
 * No, This is necessary because primitive types are not default initialized.
 *
 */
template<typename DataType>
Matrix<DataType>::Matrix() :
		number_of_rows(0), number_of_columns(0) {
}

/**
 * STP: I realized that number_of_* gets assigned again in resize. This is ok
 * but could be taken out. Is the current way as clear as possible?
 */
template<typename DataType>
Matrix<DataType>::Matrix(size_t number_of_rows, size_t number_of_columns,
		const DataType& value) :
		number_of_rows(number_of_rows), number_of_columns(number_of_columns) {
	resize(number_of_rows, number_of_columns, value);
}

template<typename DataType>
Matrix<DataType>::Matrix(const Matrix& matrix) :
		number_of_rows(matrix.number_of_rows), number_of_columns(
				matrix.number_of_columns), rows(matrix.rows) {
}

template<typename DataType>
void Matrix<DataType>::fill(const DataType& value) {
	for (size_t row = 0; row < number_of_rows; row++) {
		for (size_t column = 0; column < number_of_columns; column++)
			rows[row][column] = value;
	}
}

template<typename DataType>
void Matrix<DataType>::clear() {
	for (size_t row = 0; row < number_of_rows; ++row)
		rows[row].clear();

	number_of_rows = 0;
	number_of_columns = 0;
}

template<typename DataType>
void Matrix<DataType>::resize(size_t number_of_rows, size_t number_of_columns,
		const DataType& inVal) {
	this->number_of_rows = number_of_rows;
	this->number_of_columns = number_of_columns;
	rows.resize(number_of_rows);
	for (size_t row = 0; row < number_of_rows; row++)
		rows[row].resize(number_of_columns, inVal);
}

template<typename DataType>
inline const typename Matrix<DataType>::RowType&
Matrix<DataType>::operator[](size_t index) const {
	return rows[index];
}

template<typename DataType>
inline typename Matrix<DataType>::RowType&
Matrix<DataType>::operator[](size_t index) {
	return rows[index];
}

template<typename DataType>
inline DataType&
Matrix<DataType>::at(size_t row_index, size_t column_index) {
	if (row_index >= number_of_rows or column_index >= number_of_columns) {
		std::cerr << "ERROR! Matrix bounds exceeded" << std::endl;
		exit(1);
	}
	return rows.at(row_index).at(column_index);
}

template<typename DataType>
inline DataType&
Matrix<DataType>::at(size_t row_index, size_t column_index) const {
	if (row_index >= number_of_rows or column_index >= number_of_columns) {
		std::cerr << "ERROR! Matrix bounds exceeded" << std::endl;
		exit(-1);
	}
	return rows.at(row_index).at(column_index);
}

template<typename DataType>
inline const Matrix<DataType>&
Matrix<DataType>::operator+(const Matrix& matrix) {

	if (number_of_rows != matrix.number_of_rows) {
		std::cerr << "Number of rows are not equal" << std::endl;
		exit(-1);
	}

	if (number_of_columns == matrix.number_of_columns) {
		std::cerr << "Number of columns are not equal" << std::endl;
		exit(-1);
	}

	for (size_t row = 0; row < number_of_rows; row++) {
		rows.at(row) += matrix.rows.at(row);
	}

	return *this;
}

template<typename DataType>
inline const Matrix<DataType>&
Matrix<DataType>::operator+=(const Matrix& matrix) {
	return (*this + matrix);
}

template<typename DataType>
inline
void Matrix<DataType>::SwapRows(size_t row_1, size_t row_2) {
	if (row_1 != row_2)
		rows[row_1].swap(rows[row_2]);
}

template<typename DataType>
inline
void Matrix<DataType>::SwapColumns(size_t column_1, size_t column_2) {
	for (size_t row = 0; row < number_of_rows; row++)
		std::swap(rows[row][column_1], rows[row][column_2]);
}

template<typename DataType>
inline
void Matrix<DataType>::swap(Matrix& matrix) {
	std::swap(matrix.number_of_rows, number_of_rows);
	std::swap(matrix.number_of_columns, number_of_columns);

	rows.swap(matrix.rows);
}

template<typename DataType>
void Matrix<DataType>::symmetrize() {
	for (int row = 0; row < number_of_rows; row++) {
		for (int column = 0; column < row; column++) {
			at(row, column) = at(column, row) = (at(row, column)
					+ at(column, row)) / 2;
		}
	}
}

#endif
