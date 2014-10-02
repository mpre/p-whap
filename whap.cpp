#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <bitset>
#include <mpi.h>

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
  int my_rank, numprocs;
  int n;       //Number of fragments
  int m;       //Number of columns


  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

  vector< vector< bool > > bips; 

  if (my_rank == 0)   //MASTER PROCESSOR
    {
      ifstream ifs(argv[1], ios_base::binary);
      ifs.read((char *)&n, sizeof(int));
      ifs.read((char *)&m, sizeof(int));
      cerr << "Nrows : " << n << endl;
      cerr << "Ncols : " << m << endl;
  
      Matrix input(n, m);
      readMatrix(input, ifs);
      ifs.close();
      
      //cout << "GUARDA: " << input.toVector().size() << endl;

      // TEMPORARY: Print input matrix
      for(int row = 0; row < input.rows_num(); ++row)
        {
          for(int col = 0; col < input.cols_num(); ++col)
            {
              if(input.get(row, col) != -1)
                cerr << " " << input.get(row, col);
              else
                cerr << " -";
            }
          cerr << endl;
        }

      //Master sends to all the number of fragments and columns
      int send_nm[2];
      send_nm[0] = n;
      send_nm[1] = m;

      //Just proc 0 writes on log file in this case

      cerr << "Computing bipartitions...";

      bips = computeBipartitions(n);

      cerr << "done." << endl;
      for(int i =0; i < bips.size(); ++i)
        {
          cerr << "Bipartition #" << i << " = ";
          printBipartition(bips[i]);
          cerr << endl;
        }


      MPI_Bcast(send_nm, 2, MPI_INTEGER, 0,  MPI_COMM_WORLD);

      Matrix optimum(bips.size(), m);

      cerr << "Number of bipartitions : " << bips.size() << endl;
   
      int interval = bips.size() / (numprocs - 1); //Size of the column split
      int count = 0;
      int dest = 1;  //Rank of the receiver
      vector<int> send_opt;  //Column split


      //Compute length of the intervals

      int lengths[numprocs];
      
      //std::cerr << "LA LENGTH DI 0!!! "  << std::endl;
      for (int i = 1; i < numprocs; i++)
        {
          lengths[i] =  bips.size() / (numprocs - 1);
          if ( i <  bips.size() % (numprocs - 1))
            {
              ++lengths[i];
            }
          std::cerr << lengths[i] << " ";
        }
      std::cerr << std::endl;

      //BASE CASE
      
      MPI_Bcast(&(input.get_col(0)).front(), input.get_col(0).size(), MPI_INTEGER, 0, MPI_COMM_WORLD);
      //std::cout << my_rank << "  : INVIATO L'INPUT - CASO BASE :   " << 0 << std::endl;

      //Merge of results
      MPI_Status status;
          
      //Initializing the optimum column to the maximum value n
      for(int i = 0; i < bips.size(); i++)
        {
          optimum.set(i, 0, n);
        }


       //Receive the result from any proc
      for(dest = 1; dest < numprocs; dest++)
        {
          vector<int> result(bips.size(), 0);
          MPI_Recv(&result[0], bips.size(), MPI_INTEGER, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
          //std::cout << my_rank << "  : RICEVUTO L'OTTIMO - CASO BASE :   " << 0 << std::endl;
          
          //Update column of optimum
          for(int i = 0; i < bips.size(); i++)
            {
              if(result[i] < optimum.get(i, 0))
                {
                  optimum.set(i, 0, result[i]);
                }
            }
        }


      //ITERATIVE STEPS
      for (int col = 1; col < input.cols_num(); col++)
        {  
          MPI_Barrier(MPI_COMM_WORLD);
          count = 0;
          dest = 1;
          send_opt.clear();

          //Send the input column col to all
          MPI_Bcast(&(input.get_col(col)).front(), input.get_col(col).size(), MPI_INTEGER, 0, MPI_COMM_WORLD);
          //std::cout << my_rank << "  : INVIO LA COLONNA DI INPUT :   " << col << std::endl;

          //Send to each proc the corresponding split of optimum
          //CHECK THIS CAREFULLY
          for(int i = 0; i < bips.size(); i++)
            {
              send_opt.push_back(optimum.get(i, col - 1));
              count++;
              if (count == lengths[dest])
                {
                  MPI_Send(&send_opt.front(), send_opt.size(), MPI_INTEGER, dest, 0, MPI_COMM_WORLD);
                  //std::cout << my_rank << "  : INVIATO IL FRAMMENTO SU CUI CALCOLARE :   " << col << std::endl;
                  dest++;
                  count = 0;
                  send_opt.clear();
                }
            }

          //Merge of results
          
          //Initializing the optimum column to the maximum value n
          
          for(int i = 0; i < bips.size(); i++)
            {
              optimum.set(i, col, n);
            }
          

          //Receive the result from any proc
          for(dest = 1; dest < numprocs; dest++)
            {
              //std::cout << my_rank << "  : STO ASPETTANDO IL RISULTATO :   " << col << std::endl;
              vector<int> result(bips.size(), 0);
              MPI_Recv(&result[0], bips.size(), MPI_INTEGER, dest, 0, MPI_COMM_WORLD, &status);

              //std::cout << my_rank << "  : RICEVUTO L'OTTIMO CALCOLATO :   " << col << std::endl;

              //Update column of optimum
              for(int i = 0; i < bips.size(); i++)
                {
                  if(result[i] < optimum.get(i, col))
                    {
                      optimum.set(i, col, result[i]);
                    }
                }
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

    } else { //NOT MASTER PROCESSORS

    MPI_Status status;
    

    //Each proc receive the number of fragments and columns
    int result[2];
    MPI_Bcast(result, 2, MPI_INTEGER, 0,  MPI_COMM_WORLD);
    
    n = result[0];
    m = result[1];

    bips = computeBipartitions(n);

    //Compute length of the intervals

    int lengths[numprocs];
    int starts[numprocs];
      
    for (int proc = 1; proc < numprocs; proc++)
      {
        lengths[proc] =  bips.size() / (numprocs - 1);
        if ( proc <  bips.size() % (numprocs - 1))
          {
            ++lengths[proc];
          }
        for (int k = 1; k < proc; k++)
          {
            starts[proc] += lengths[k];
          }
      }
    
    //Resulting vector
    vector<int> optimum_col_new;

    //BASE CASE
        
    //Receive the input column
    vector<int> input(n, 0);
    MPI_Bcast(&input[0], n, MPI_INTEGER, 0, MPI_COMM_WORLD);
    //std::cout << my_rank << "  : RICEVUTO L'INPUT - CASO BASE :   " << 0 << std::endl;
        
    int delta = 0;   // Local contribution to opt solution
    vector< int > act_pos = computeActivePositions(input);

    optimum_col_new.clear();
            
    //Initializing the resulting vector with the maximum values n
    for(int i = 0; i < bips.size(); i++)
      {
        optimum_col_new.push_back(n);
      }

    //Computing Case Base values
    for(int row = 0; row < bips.size(); row++)
      {
        delta = computeDelta(input, act_pos, bips[row]);
        optimum_col_new[row] = delta;
      }
    //Remembering the input of the previous step
    vector<int> input_prec = input;
            
    //Send the result
    MPI_Send(&optimum_col_new.front(), optimum_col_new.size(), MPI_INTEGER, 0, 0, MPI_COMM_WORLD);
    //std::cout << my_rank << "  : INVIATO L'OTTIMO - CASO BASE :   " << 0 << std::endl;

    //ITERATIVE STEPS

    for(int col = 1; col < m; col++)
      {
        MPI_Barrier(MPI_COMM_WORLD);
        
        // Receive the column of the input
        vector<int> input(n, 0);
        MPI_Bcast(&input[0], n, MPI_INTEGER, 0, MPI_COMM_WORLD);
        //std::cout << my_rank << "  : RICEVUTO INPUT NUOVO :   " << col << std::endl;

        delta = 0;   // Local contribution to opt solution
        int minimum = 0; // Min value of according bipartition in col-1
        vector< int > act_pos = computeActivePositions(input);

        optimum_col_new.clear();
            
        //Initializing the resulting vector with the maximum values n
        
        for(int i = 0; i < bips.size(); i++)
          {
            optimum_col_new.push_back(n);
          }
        
        //Receive the fragment of the optimum in the previous step
        //std::cout << my_rank << "  : STO ASPETTANDO IL FRAMMENTO SU CUI LAVORARE :   " << col << std::endl;
        vector<int> optimum_fragment(lengths[my_rank], 0);
        MPI_Recv(&optimum_fragment[0], lengths[my_rank], MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //std::cout << my_rank << "  : RICEVUTO FRAMMENTO SU CUI LAVORARE :   " << col << std::endl;

        vector<int> optimum_prec;
        

        // DA SISTEMAREEEEEEEEEEEEEEEEEEEEEEEEE!!!!!
        int index = 0;
        
        for(int i = 0; i < bips.size(); i++)
          {
            if(i < starts[my_rank] || i >= starts[my_rank] + lengths[my_rank])
              {
                optimum_prec.push_back(n);
              } else {
              optimum_prec.push_back(optimum_fragment[index]);
              index++;
            }            
          }
            
        for(int row = 0; row < bips.size(); row++)
          {
            delta = computeDelta(input, act_pos, bips[row]);
            minimum = computeMinimum(input_prec, optimum_prec,
                                     act_pos, bips[row], bips);
            optimum_col_new[row] = delta + minimum;
          }
            
        //Remembering the input of the previous step
        input_prec = input;
            
        //Send the result
        MPI_Send(&optimum_col_new.front(), optimum_col_new.size(), MPI_INTEGER, 0, 0, MPI_COMM_WORLD);

        //std::cout << my_rank << "  : INVIATO IL MIO OTTIMO :   " << col << std::endl;
      } 
  }

  MPI_Finalize();
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
        //cerr << b_index << " ";
        if(opt_col[b_index] < minimum)
          minimum = opt_col[b_index];
      }

  return minimum;
}

int computeDelta(const vector< int >& frag_col, const vector< int >& act_pos,
                 const vector< bool >& cbip)
{
  enum {OO, OI, IO, II};
  vector< int > solutions(II+1, 0);

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
