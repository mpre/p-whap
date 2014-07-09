#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <algorithm>

#include "Matrix.hpp"
#include "Bipartition.hpp"

using std::vector;

//Definition of the functions

/*
  Compute  the  intersection  between  act_pos_1 and  act_pos_2  and  build  the
  corrisponding subsequences on bip1 and bip2.

  Return true if they are equal or complementary.
  Return false otherwise.

  Additional info:
  This function check whether the active  positions (its number equal to len) on
  the bipartition bip1 and  the active positions (its number equal  to len 2) on
  the bipartition bip2 are bipartited in the same way.
*/
bool accordance(const vector<bool>& bip1, const vector<int>& act_pos_1,
                const vector<bool>& bip2, const vector<int>& act_pos_2);

/*
  Return the suboptimal value for cbip.

  Additional info:
  Given the  column j, with  the corresponding  active positions (its  number is
  equal  to len)  activePoitionsJ and  given the  bipartition bipJ,  we want  to
  compute minimum  that is  the minimum  value in  a column  j-1 of  the optimum
  matrix corresponding to a bipartition bip that is in accordance with bipJ.
*/
int computeMinimum(const vector< int >& frag_col, const vector< int >& opt_col,
                   const vector< int >& active_pos, const vector< bool >& cbip,
                   const vector< vector< bool > >& bip_set);

/*
  Return the  minimum number of corrections  to turn frag_col into  a homozygous
  column or a column in accordance to cbip.
  
  Additional info:
  Given a  column col,  its active positions  active with a  length of  len, and
  given a bipartition bip  we want to compute the minimum  number of changes for
  the 4  possible combinations (both parts  equal to 0, or  one 0 and one  1, or
  both equal to 1) corresponding to the given bipartition.
*/
int computeDelta(const vector< int >& frag_col, const vector< int >& act_pos,
                 const vector< bool >& cbip);

/*
  Return the active position in frag_col.

  Additional info:
  Given  a  column j,  this  function  compute  the array  activePositions  that
  contains all the positions  (i.e. fragments) covered by the column  j and l is
  the length of this computed array.
*/
vector< int > computeActivePositions(const vector< int >& frag_col);

/* Return all the possible bipartitions (max frag_num is 2^64)*/
void computeBipartitions(vector< vector<bool> >& bips, const int frags_num);

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
  int len;
  
  Matrix input(n, m);
  Matrix optimum(b, m);
  
  vector< vector< bool > > bips;

  computeBipartitions(bips, n);

  // Inserire caso base
  // Per opt calcolare solo delta e riempire la prima colonna
  for(int col = 1; col < m; col++)
    {
      int delta = 0;   // Local contribution to opt solution
      int minimum = 0; // Min value of according bipartition in col-1
      vector< int > act_pos = computeActivePositions(input.get_col(col));

      for(int row = 0; row < b; row++)
        {
          delta = computeDelta(input.get_col(col), act_pos, bips[row]);
          //          computeDelta(col, activePositions, len, bips[row], delta, input);

          minimum = computeMinimum(input.get_col(col-1), optimum.get_col(col-1),
                                   act_pos, bips[row], bips);
          //          computeMinimum(col, activePositions, len, bips[row], minimum, optimum, input, bips);
          optimum.set(row, col, delta + minimum);
        }
    }

  return 0;
}

bool accordance(const vector<bool>& bip1, const vector<int>& act_pos_1,
                const vector<bool>& bip2, const vector<int>& act_pos_2)
{
  // Compute shared_positions
  vector< int > shared_pos(std::max(act_pos_1.size(), act_pos_2.size()));
  vector< int >::iterator it = std::set_intersection( act_pos_1.begin(), act_pos_1.end(),
                                                      act_pos_2.begin(), act_pos_2.end(),
                                                      shared_pos.begin() );
  shared_pos.resize(it - shared_pos.begin());

  bool complemented = true;
  bool equal = true;
  for(int i =0; i < shared_pos.size(); ++i)
    {
      if(bip1[shared_pos[i]] == bip2[shared_pos[i]])
        complemented = false;
      else
        equal = false;
    }

  return (equal || complemented);
}

int computeMinimum(const vector< int >& frag_col, const vector< int >& opt_col,
                   const vector< int >& active_pos, const vector< bool >& cbip,
                   const vector< vector< bool > >& bip_set)
{
  vector< int > prev_act_pos = computeActivePositions(frag_col);

  int minimum = frag_col.size();

  for(int b_index = 0; b_index < bip_set.size(); b_index++)
    if(accordance(cbip, active_pos, bip_set[b_index], prev_act_pos))
      if(opt_col[b_index] < minimum)
        minimum = opt_col[b_index];

  return minimum;
}

int computeDelta(const vector< int >& frag_col, const vector< int >& act_pos,
                 const vector< bool >& cbip)
{
  enum {OO, OI, IO, II};
  vector< int > solutions(II, 0);

  for(int p = 0; p < act_pos.size(); ++p)
    {
      if(cbip[act_pos[p]] == false)
        {
          // If the element active[i] is in the part 0 in the bipartition bip
          if(frag_col[act_pos[p]] == 1)
            {
              ++solutions[OO];
              ++solutions[OI];
            }
          else
            {
              ++solutions[IO];
              ++solutions[II];
            }
        }
      else
        {
          // If the element active[i] is in the part 1 in the bipartition bip
          if (frag_col[act_pos[p]] == 1) {
            ++solutions[OO];
            ++solutions[IO];
          } else {
            ++solutions[OI];
            ++solutions[II];
          }
        }
    }

  return std::min( std::min(solutions[OO], solutions[OI]),
                   std::min(solutions[IO], solutions[II]));
}

vector< int > computeActivePositions(const vector< int >& frag_col)
{
  vector< int > act_pos;
  for(int i =0; i == frag_col.size(); ++i)
    if(frag_col[ i ] != -1)
      act_pos.push_back(i);
}

void computeBipartitions(vector< vector<bool> > &bips, int frags_num)
{
  for(int i =0; i < pow(2, frags_num); ++i)
    {
      vector< bool > cbip(frags_num, false);
      for(int j =0; j < frags_num; ++j)
        {
          cbip[j] = (i & (1 << j));
        }
      bips.push_back(cbip);
    }
}

