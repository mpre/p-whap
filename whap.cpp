#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <bitset>

#include "Matrix.hpp"
#include "Bipartition.hpp"

using namespace std;

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
vector< vector< bool > > computeBipartitions(const int frags_num);

// Support functions

/* Read matrix from binary ifstream */
void readMatrix(Matrix& input, ifstream& ifs);

/* Print a Matrix instance */
void printMatrix(Matrix& in, vector< vector< bool > >& bips);

/* Print bipartition */
void printBipartition(vector< bool >& bip);

//MAIN

int main(int argc, char** argv)
{	
  int n;       //Number of fragments
  int m;       //Number of columns

  ifstream ifs(argv[1], ios_base::binary);
  ifs.read((char *)&n, sizeof(int));
  ifs.read((char *)&m, sizeof(int));
  cerr << "Nrows : " << n << endl;
  cerr << "Ncols : " << m << endl;
  
  Matrix input(n, m);
  readMatrix(input, ifs);
  ifs.close();

  // TEMPORARY: Print input matrix
  for(int row =0; row < input.rows_num(); ++row)
    {
      for(int col =0; col < input.cols_num(); ++col)
        {
          if(input.get(row, col) != -1)
            cerr << " " << input.get(row, col);
          else
            cerr << " -";
        }
      cerr << endl;
    }
  cerr << "Computing bipartitions...";
  vector< vector< bool > > bips = computeBipartitions(n);
  cerr << "done." << endl;

  for(int i =0; i < bips.size(); ++i)
    {
      cerr << "Bipartition #" << i << " = ";
      printBipartition(bips[i]);
      cerr << endl;
    }

  Matrix optimum(bips.size(), m);

  cerr << "Number of bipartitions : " << bips.size() << endl;
  
  for(int i =0; i < bips.size(); ++i)
    {
      vector< int > act_pos = computeActivePositions(input.get_col(0));
      int opt_val = computeDelta(input.get_col(0), act_pos, bips[i]);
      optimum.set(i, 0, opt_val);
    }

  // Inserire caso base
  // Per opt calcolare solo delta e riempire la prima colonna
  for(int col = 1; col < input.cols_num(); col++)
    {
      cerr << "Step " << col << "/" << input.cols_num() << endl;
      int delta = 0;   // Local contribution to opt solution
      int minimum = 0; // Min value of according bipartition in col-1
      vector< int > act_pos = computeActivePositions(input.get_col(col));

      for(int row = 0; row < optimum.rows_num(); row++)
        {
          delta = computeDelta(input.get_col(col), act_pos, bips[row]);
          cerr << "Bip #" << row << " in accordance with: ";
          minimum = computeMinimum(input.get_col(col-1), optimum.get_col(col-1),
                                   act_pos, bips[row], bips);
          cerr << endl;
          optimum.set(row, col, delta + minimum);
        }
    }

  printMatrix(optimum, bips);

  vector< int > last_col = optimum.get_col(m -1);
  int best_index = min_element(last_col.begin(), last_col.end()) - last_col.begin();
  cerr << "Optimum is " << *min_element(last_col.begin(), last_col.end())
       << " found at position "
       << best_index << endl;

  printBipartition(bips[best_index]);
  cerr << endl;

  return 0;
}

bool accordance(const vector<bool>& bip1, const vector<int>& act_pos_1,
                const vector<bool>& bip2, const vector<int>& act_pos_2)
{
  // Compute shared_positions
  vector< int > shared_pos(max(act_pos_1.size(), act_pos_2.size()));
  vector< int >::iterator it = set_intersection( act_pos_1.begin(), act_pos_1.end(),
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
      {
        cerr << b_index << " ";
        if(opt_col[b_index] < minimum)
          minimum = opt_col[b_index];
      }

  return minimum;
}

int computeDelta(const vector< int >& frag_col, const vector< int >& act_pos,
                 const vector< bool >& cbip)
{
  enum {OO, OI, IO, II};
  vector< int > solutions(II, 0);

  for(int p = 0; p < act_pos.size(); ++p)
    {
      if(cbip[act_pos[p]] == false) // Posizione della bipartizione attiva a 0
        {
          // If the element active[i] is in the part 0 in the bipartition bip
          if(frag_col[act_pos[p]] == 1)
            {
              solutions[OO] += 1;
              solutions[OI] += 1;
            }
          else
            {
              solutions[IO] += 1;
              solutions[II] += 1;
            }
        }
      else
        {
          // If the element active[i] is in the part 1 in the bipartition bip
          if (frag_col[act_pos[p]] == 1) {
            solutions[OO] += 1;
            solutions[IO] += 1;
          } else {
            solutions[OI] += 1;
            solutions[II] += 1;
          }
        }
    }

  return *min_element(solutions.begin(), solutions.end());
}

vector< int > computeActivePositions(const vector< int >& frag_col)
{
  vector< int > act_pos;
  for(int i =0; i < (int)frag_col.size(); ++i)
    if(frag_col[ i ] != -1)
      act_pos.push_back(i);

  return act_pos;
}

vector< vector< bool > > computeBipartitions(const int frags_num)
{
  vector< vector< bool > > bips;
  for(int i =0; i < pow(2, frags_num); ++i)
    {
      vector< bool > cbip(frags_num, false);
      for(int j =0; j < frags_num; ++j)
        {
          cbip[j] = (i & (1 << j));
        }
      bips.push_back(cbip);
    }
  return bips;
}

void printMatrix(Matrix& in, vector< vector< bool > >& bips)
{
  for(int row = 0; row < in.rows_num(); ++row)
    {
      printBipartition(bips[row]);
      cerr << ": ";
      for(int col =0; col < in.cols_num(); ++col)
        {
          cerr << in.get(row, col) << " ";
        }
      cerr << endl;
    }
}

void readMatrix(Matrix& input, ifstream& ifs)
{
  for(int row =0; row < input.rows_num(); ++row)
    for(int col =0; col < input.cols_num(); ++col)
      {
        short t;
        ifs.read((char *)&t, sizeof(short));
        input.set(row, col, t);
      }
}

void printBipartition(vector< bool >& bip)
{
  for(vector<bool>::iterator it = bip.begin(); it != bip.end(); ++it)
    cerr << *it;
}
