#include <Bipartition.hpp>

Bipartition::Bipartition(int *new_bip, int length)
{
  bip = new int[length];
  for (int i = 0; i < length; i++) 
    {
      bip[i] = new_bip[i];
    }
  len = length;
}

Bipartition::~Bipartition()
{
  delete[] bip;
}
  
int Bipartition::at(int pos)
{
  return bip[pos];
}


