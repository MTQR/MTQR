#include "mtqr.h"

int main(int argc, char** argv)
{
  // Initialise the values for the input parameters of the silent mode
  double lambda_min = static_cast<double>(-1/E);
  double lambda_max = static_cast<double>(1.0/2.0);

  // Instantiate MTQR's primary module in silent mode
  std::vector<std::vector<double>> array = mtqr(lambda_min, lambda_max);

  return 0;
}