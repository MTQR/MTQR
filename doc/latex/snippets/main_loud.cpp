#include "mtqr.h"

int main(int argc, char** argv)
{
  std::vector<float128> coeff_sequence = {E, 5.0, -1.0, 1.0, 10.0};
  std::vector<float128> muntz_sequence = {E + 0.25, -PI/4, -0.5, 0, 2};
  mtqr(muntz_sequence, coeff_sequence);

  return 0;
}