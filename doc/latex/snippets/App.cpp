#include "Quasimont.h"

int main()
{
    std::vector <float128> coeff_sequence = {1, 1, 1, 1, 1};
  std::vector <float128> muntz_sequence = {0,
                       1,
                       boost::math:: constants ::pi<float128 >()/2,
                       boost::math:: constants ::pi<float128 >()/3,
                       boost::math:: constants ::pi<float128 >() /6};

    quasimont(muntz_sequence, coeff_sequence);

    // Da inserire qui il loading dei pesi e nodi trasformati da quasimont attraverso
    // i files output/Nodes.txt e Weights.txt

  return 0;
}
