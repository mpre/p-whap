#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>

using std::string;

static const int NUM_PARAMS =6;
static const short emptypos   =-1;

void print_help()
{
  std::cerr << "Usage: inputgen ";
  std::cerr << "<frags>  ";
  std::cerr << "<flen> ";
  std::cerr << "<slen> ";
  std::cerr << "<outputfilename> ";
  std::cerr << std::endl;
  std::cerr << std::endl << "Options:" << std::endl;
  std::cerr << "\t  <frags>              # numbers of fragments" << std::endl;
  std::cerr << "\t  <flen>               # length of the fragments" << std::endl;
  std::cerr << "\t  <slen>               # length of the sequence" << std::endl;
  std::cerr << "\t  <outputfilename>     # filename" << std::endl;
  std::cerr << std::endl;
  std::cerr << std::endl << "Notes:" << std::endl;
  std::cerr << "\t fragments will have length = len/frags" << std::endl;
}

int main(int argc, char** argv)
{
  if(argc != NUM_PARAMS)
    {
      print_help();
      return -1;
    }

  int slen      = atoi(argv[1]);
  int flen      = atoi(argv[2]);
  int blocknum  = atoi(argv[3]);
  int blockwidth= atoi(argv[4]);

  string outputfilename(argv[5]);
  
  int step = std::floor((float) slen / blocknum);
  
  std::cerr << "Length of sequence  : " << slen << std::endl;
  std::cerr << "Length of framents  : " << flen << std::endl; 
  std::cerr << "Number of blocks    : " << blocknum << std::endl;
  std::cerr << "Block width         : " << blockwidth << std::endl;
  std::cout << "Step                : " << step << std::endl;
  std::cerr << "Outputfilename      : " << outputfilename << std::endl << std::endl;

  if(step > flen || flen*blocknum < slen)
    {
      std::cerr << "[ERROR] Step too wide, aborting." << std::endl;
      return -1;
    }

  std::cout << "Fragments will cover..." << std::endl;
  for(int i=0; i < blocknum; i+=step)
    {
      std::cout << "[ " << std::min(i, slen - flen) << ", "
                << std::min(i + flen, slen) << " ]" << std::endl;
    }

  std::ofstream output(outputfilename.c_str(), std::ios_base::out | std::ios_base::binary);
  // Write number of rows
  int nrows = blocknum * blockwidth;
  output.write((char *)&nrows, sizeof(int));
  // Write number of columns
  output.write((char *)&slen, sizeof(int));

  for(int i=0; i < blocknum; ++i)
    {
      for(int j =0; j < blockwidth; ++j)
        {
          short outputc = (j % 2) ? 1 : 0;
          std::cout << "Frags (" << i << ") #" << j << ":";
          for(int k = 0; k <slen; ++k)
            {
              if(k >= std::min(i * step, slen - flen) && k < std::min(i*step + flen, slen))
                {
                  std::cout << " " << outputc;
                  output.write((char *)&outputc, sizeof(short));
                }
              else
                {
                  std::cout << " -";
                  output.write((char *)&emptypos, sizeof(short));
                }
            }
          std::cout << std::endl;
        }
    }

  output.close();

  return 0;
}
