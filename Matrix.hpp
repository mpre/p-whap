#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

using std::vector;

class Matrix 
{

private: 

  int rows;
  int columns;
  
  int* matrix;

public:
  
  Matrix();

  Matrix(int rows_in, int columns_in);
  
  Matrix(int rows_in, int columns_in, const vector<int>& vect);

  ~Matrix();

  void set(int row, int col, int value);
  
  int get(int row, int col);

  vector< int > get_col(int index);

  int cols_num() const;

  int rows_num() const;
  
  vector<int> toVector();
};

#endif
