#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//Data type

/* 
   This is the type of data necessary to represent a bi-partization of the fragments. 
   It has to represent an array with length equal to the number of fragments
   where each position i of this array corresponds to fragment i and its value is equal
   to 0 or to 1, depending on which part fragment i has been placed.

   Notice that two bi-partition are equivalent if they are equal or complementary. 
*/
typedef int Bipartition;


/*
  This is the type necessary to represent a fragment. Its length is equal to the number of
  SNP positions in the input. Each position i is equal to a simbol "-" if the fragment does not
  cover the corresponding position i, or it is equal to 0/1 depending on value in the input.
  A fragment corresponds to a sequence of 0's and 1's, preceded and followed by a certain number
  of symbols "-". Between each pair of 0's and 1's in a fragment there cannot be any symbol "-".
*/
typedef int Fragment;


/*
  This is the type necessary to representa a matrix (n x m) with a certain number of rows n 
  and a certain number of columns m
*/
typedef int Matrix;



//Definition of the functions


bool accordance(Bipartition bip1, int* active1, int len, Bipartition bip2, int* active2, int len2);
void computeMinimum(int j, int* activePositionsJ, int len, Bipartition bipJ, int &minimum);
void computeDelta(int col, int* active, int len, Bipartition bip, int &delta);
void computeActivePositions(int col, int &l, int *activePositions);
void computeBipartitions(Bipartition *bips);
bool isIn(int n, int *array, int len);
int at(Bipartition bip, int position);
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
Matrix input;


/*
  This is the matrix (b x m) that will contain the suboptimal structure, such that b is the number
  of all the possible bipartitions (notice that b = 2^(n) where n is still the number of fragments)
  and m is the number of positions. For each bipartition i and each position j, then 
  optimum[i, j] is equal to the optimal solution for the columns 1..j using the bipartition i on
  the last column j.
*/
int b;       //Number of bipartitions
Matrix optimum;


/*
  Set of all the possible bipartitions
*/

Bipartition *bips;




//MAIN



int main(int argc, char** argv)
{	
  int delta;
  int minimum;
  int* activePositions;
  int len;

  
  computeBipartitions(bips);


  for (int col = 0; col < m; col++)
    {
      computeActivePositions(col, len, activePositions);
      delta = 0;
      minimum = 0;

      for(int i = 0; i < b; i++)
        {
          computeDelta(col, activePositions, len, bips[i], delta);
          computeMinimum(col, activePositions, len, bips[i], minimum);

          setMatrix(optimum, i, col, delta + minimum);
        }

    }

}


/*
  This function check whether the active positions (its number equal to len) on the bipartition bip1
  and the active positions (its number equal to len 2) on the bipartition bip2 are bipartited in the
  same way.
*/

bool accordance(Bipartition bip1, int* active1, int len, Bipartition bip2, int* active2, int len2)
{
  int pos = 0;
  bool flagEqual = true;
  bool flagComplementary = true;

  for(int i = 0; i < len; i++)
    {
      if (isIn(active1[i], active2, len2))
        {
          pos = active1[i];
          if (at(bip1, pos) == at(bip2, pos))
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


void computeMinimum(int j, int* activePositionsJ, int len, Bipartition bipJ, int &minimum)
{
  int len1;
  int* activePositionsJ1;
  int tempMinimum = 0;

  computeActivePositions(j-1, len1, activePositionsJ1);

  for (int i = 0; i < b; i++)
    {
      if (accordance(bipJ, activePositionsJ, len, bips[i], activePositionsJ1, len1))
        {
          if(getMatrix(optimum, i, j - 1) < tempMinimum)
            {
              tempMinimum = getMatrix(optimum, i, j - 1);
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

void computeDelta(int col, int* active, int len, Bipartition bip, int &delta)
{
  //Value for the combination when both the parts are equal to 0
  int sol1 = 0;
  //Value for the combination when part 0 is equal to 0 and part 1 equal to 1
  int sol2 = 0;
  //Value for the combination when part 0 is equal to 1 and part 1 equal to 0
  int sol3 = 0;
  //Value for the combination when both the parts are equal to 1
  int sol4 = 0;

  int i = 0;

  for(i = 0; i < len; i++)
    {
      //If the element active[i] is in the part 0 in the bipartition bip
      if (at(bip, active[i]) == 0)
        {
          if (getMatrix(input, active[i], col) == 1) {
            sol1++;
            sol2++;
          } else {
            sol3++;
            sol4++;
          }
          //If the element active[i] is in the part 1 in the bipartition bip
        } else {
        if (getMatrix(input, active[i], col) == 1) {
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


void computeActivePositions(int col, int &l, int *activePositions) {
  int len = 0;
  int temp;

  for (int i = 0; i < m; i++) {
    temp = getMatrix(input, i, col);
    if(temp == 0 || temp == 1) {
      len++;
    }
  }

  int tempActive[len];
  int iter = 0;
  
  for (int i = 0; i < m; i++) {
    temp = getMatrix(input, i, col);
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


void computeBipartitions(Bipartition *bips){}


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


/*
  Function that returns the element of the bipartition bip in the position "position"
*/

int at(Bipartition bip, int position) {}

/*
  Function that sets the value "value" at position (i,j) of the Matrix matrix.
*/

void setMatrix(Matrix matrix, int i, int j, int value) {}


/*
  Function that returns the value of the Matrix matrix in the position (i,j)
*/
int getMatrix(Matrix matrix, int i, int j) {}
