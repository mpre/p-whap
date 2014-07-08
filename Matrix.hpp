#ifndef MATRIX_H
#define MATRIX_H

class Matrix 
{

private: 

  int rows;
  int columns;
  
  int* matrix;

public:
  
  Matrix(int rows_in, int columns_in);

  ~Matrix();

  void setMatrix(int i, int j, int value);
  
  int getMatrix(int i, int j);
};

#endif
