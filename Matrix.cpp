#include <Matrix.hpp>  

Matrix::Matrix()
{
  rows = 0;
  columns = 0;
}

Matrix::Matrix(int rows_in, int columns_in)
{
  rows = rows_in;
  columns = columns_in;
  matrix.resize(rows * columns);
}

Matrix::Matrix(int rows_in, int columns_in, const vector<int>& vect)
{
  rows = rows_in;
  columns = columns_in;
  matrix.resize(rows * columns);
  for(int i = 0; i < vect.size(); i++)
    {
      matrix[i] = vect[i];
    }
}

/**
Matrix::~Matrix()
{
  delete[] matrix;
}
**/

void Matrix::init(int rows_in, int columns_in)
{
  rows = rows_in;
  columns = columns_in;
  matrix.resize(rows * columns);
}

void Matrix::init(int rows_in, int columns_in, const vector<int>& vect)
{
  rows = rows_in;
  columns = columns_in;
  matrix.resize(rows * columns);
  for(int i = 0; i < vect.size(); i++)
    {
      matrix[i] = vect[i];
    }
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

int Matrix::cols_num() const
{
  return columns;
}

int Matrix::rows_num() const
{
  return rows;
}


vector<int> Matrix::toVector()
{
  vector<int> res;
  for(int i = 0; i < rows * columns; i++)
    {
      res.push_back(matrix[i]);
    }
  return res;
}
