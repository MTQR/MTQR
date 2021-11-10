#include "Quasimont.h"

int main(int argc, char** argv)
{
  std::vector<double> interval = {0, 1};
  std::vector<float128> coefficients_sequence = {1, 1, 1, 1, 1};
  std::vector<float128> muntz_sequence = {0,
                                          1,
                                          boost::math::constants::pi<float128>()/2,
                                          boost::math::constants::pi<float128>()/3,
                                          boost::math::constants::pi<float128>()/6};
  
  quasimont(muntz_sequence, coefficients_sequence, interval, plots);
  return 0;
}