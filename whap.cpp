#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Matrix.hpp>
#include <Bipartition.hpp>
#include <vector>

using std::vector;

//Definition of the functions

bool accordance(vector<int> &bip1, int *active1, int len, vector<int> &bip2, int *active2, int len2);
void computeMinimum(int j, int* activePositionsJ, int len, vector<int> &bipJ, int &minimum, Matrix &optimum, Matrix &input, vector<vector<int> > &bips);
void computeDelta(int col, int* active, int len, vector<int> &bip, int &delta, Matrix &input);
void computeActivePositions(int col, int &l, int *activePositions, Matrix &input);
void computeBipartitions(vector<vector<int> > &bips);
bool isIn(int n, int *array, int len);
void setMatrix(Matrix matrix, int i, int j, int value);
int getMatrix(Matrix matrix, int i, int j);

//Global data

/*
  This (n x m) is the input matrix, such that it contains n fragments that cover m SNP positions
  (i.e. each position is covered by at least one fragment). Each row is a fragment and each column
  is a SNP position.
*/
int n;       //Number of fragments
int m;       //Number of columns

/*
  This is the matrix (b x m) that will contain the suboptimal structure, such that b is the number
  of all the possible bipartitions (notice that b = 2^(n) where n is still the number of fragments)
  and m is the number of positions. For each bipartition i and each position j, then 
  optimum[i, j] is equal to the optimal solution for the columns 1..j using the bipartition i on
  the last column j.
*/
int b;       //Number of bipartitions

//MAIN

int main(int argc, char** argv)
{	
  int delta;
  int minimum;
  int* activePositions;
  int len;
  
  Matrix input(n, m);
  Matrix optimum(b, m);
  
  vector<vector<int> > bips;

  computeBipartitions(bips);

  for (int col = 0; col < m; col++)
    {
      computeActivePositions(col, len, activePositions, input);
      delta = 0;
      minimum = 0;

      for(int i = 0; i < b; i++)
        {
          computeDelta(col, activePositions, len, bips.at(i), delta, input);
          computeMinimum(col, activePositions, len, bips.at(i), minimum, optimum, input, bips);

          optimum.setMatrix(i, col, delta + minimum);
        }
    }
  delete activePositions;

  return 0;
}

/*
  This function check whether the active positions (its number equal to len) on the bipartition bip1
  and the active positions (its number equal to len 2) on the bipartition bip2 are bipartited in the
  same way.
*/
bool accordance(vector<int> &bip1, int *active1, int len, vector<int> &bip2, int *active2, int len2)
{
  int pos = 0;
  bool flagEqual = true;
  bool flagComplementary = true;

  for(int i = 0; i < len; i++)
    {
      if (isIn(active1[i], active2, len2))
        {
          pos = active1[i];
          if (bip1[pos] == bip2[pos])
            {
              flagComplementary = false;
            } else {
            flagEqual = false;
          }
        }
    }
  if (flagEqual || flagComplementary)
    {
      return true;
    } else {
    return false;
  }
}

/*
  Given the column j, with the corresponding active positions (its number is equal to len) activePoitionsJ
  and given the bipartition bipJ, we want to compute minimum that is the minimum value in a column j-1 of the
  optimum matrix corresponding to a bipartition bip that is in accordance with bipJ.
*/
void computeMinimum(int j, int* activePositionsJ, int len, vector<int> &bipJ, int &minimum, Matrix &optimum, Matrix &input, vector<vector<int> > &bips)
{
  int len1;
  int* activePositionsJ1;
  int tempMinimum = 0;

  computeActivePositions(j-1, len1, activePositionsJ1, input);

  for (int i = 0; i < b; i++)
    {
      if (accordance(bipJ, activePositionsJ, len, bips[i], activePositionsJ1, len1))
        {
          if(optimum.getMatrix(i, j - 1) < tempMinimum)
            {
              tempMinimum = optimum.getMatrix(i, j - 1);
            }
        }
    }
  minimum = tempMinimum;
}

/*
  Given a column col, its active positions active with a length of len, and given a bipartition
  bip we want to compute the minimum number of changes for the 4 possible combinations (both parts
  equal to 0, or one 0 and one 1, or both equal to 1) corresponding to the given bipartition.
*/
void computeDelta(int col, int* active, int len, vector<int> &bip, int &delta, Matrix &input)
{
  //Value for the combination when both the parts are equal to 0
  int sol1 = 0;
  //Value for the combination when part 0 is equal to 0 and part 1 equal to 1
  int sol2 = 0;
  //Value for the combination when part 0 is equal to 1 and part 1 equal to 0
  int sol3 = 0;
  //Value for the combination when both the parts are equal to 1
  int sol4 = 0;

  for(int i = 0; i < len; i++)
    {
      //If the element active[i] is in the part 0 in the bipartition bip
      if (bip.at(active[i]) == 0)
        {
          if (input.getMatrix(active[i], col) == 1) {
            sol1++;
            sol2++;
          } else {
            sol3++;
            sol4++;
          }
          //If the element active[i] is in the part 1 in the bipartition bip
        } else {
        if (input.getMatrix(active[i], col) == 1) {
          sol1++;
          sol3++;
        } else {
          sol2++;
          sol4++;
        }
      }
    }
  
  if (sol1 < sol2) {
    if (sol1 < sol3) {
      if (sol1 < sol4) {
        delta = sol1;
      } else {
        delta = sol4;
      }
    } else {
      if (sol3 < sol4) {
        delta = sol3;
      } else {
        delta = sol4;
      }
    }
  } else {
    if (sol2 < sol3) {
      if (sol2 < sol4) {
        delta = sol2;
      } else {
        delta = sol4;
      }
    } else {
      if (sol3 < sol4) {
        delta = sol3;
      } else {
        delta = sol4;
      }
    }
  }
}

/*
  Given a column j, this function compute the array activePositions that contains
  all the positions (i.e. fragments) covered by the column j and l is the length
  of this computed array.
*/
void computeActivePositions(int col, int &l, int *activePositions, Matrix &input) {
  int len = 0;
  int temp;

  for (int i = 0; i < m; i++) {
    temp = input.getMatrix(i, col);
    if(temp == 0 || temp == 1) {
      len++;
    }
  }

  int tempActive[len];
  int iter = 0;
  
  for (int i = 0; i < m; i++) {
    temp = input.getMatrix(i, col);
    if(temp == 0 || temp == 1) {
      tempActive[iter] = i;
      iter++;
    }
  }
  
  l = len;
  activePositions = &tempActive[len];
}

/*
  The function that computes the set of all the possible bipartitions of all the fragments.
*/
void computeBipartitions(vector<vector<int> > &bips)
{
  vector<int> temp;

  for(int j = 0; j < n; j++)
  {
    temp.push_back(0);
  }

  for (int i = 0; i < pow(2, n); i++)
  {
    for (int j = 0; j < n; j++)
    {
      if (((1 << j) & i) == 0)
      {
        temp[j] = 0;
      } else {
        temp[j] = 1;
      }
    }
    bips.push_back(temp);
  }
}

/*
  Function that returns true if and only if the array "array" with a length equal to len
  contains the value n.
*/
bool isIn(int n, int *array, int len)
{
  for(int i = 0; i < len; i++)
    {
      if (n == array[len])
        {
          return true;
        }
    }
  return false;
}

