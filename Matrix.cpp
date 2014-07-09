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

void Matrix::set(int row, int col, int value)
{
  matrix[col * rows + row] = value;
}
  
int Matrix::get(int row, int col)
{
  return matrix[col * rows + row];
}

vector< int > Matrix::get_col(int index)
{
  vector< int > currentcol;
  for(int p = index*rows; p < (index + 1)*rows; ++p)
    currentcol.push_back(matrix[p]);
  return currentcol;
}
