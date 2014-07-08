#include <Matrix.hpp>  

Matrix::Matrix(int rows_in, int columns_in)
  {
    rows = rows_in;
    columns = columns_in;
    matrix = new int[rows * columns];
  }

Matrix::~Matrix()
  {
    delete[] matrix;
  }

void Matrix::setMatrix(int i, int j, int value)
  {
    matrix[i * rows + j] = value;
  }
  
int Matrix::getMatrix(int i, int j)
  {
    return matrix[i * rows + j];
  }
