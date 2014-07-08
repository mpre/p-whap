#ifndef BIPARTITION_H
#define BIPARTITION_H

class Bipartition 
{

private: 

  int len;
  
  int* bip;

public:
  
  Bipartition(int *new_bip, int length);

  ~Bipartition();
  
  int at(int pos);
};

#endif
