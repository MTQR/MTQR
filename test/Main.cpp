//---------------------------------------------------------------------------------------
// File:      test/Main.cpp
//
// Library:   QUASIMONT-QUAdrature of SIngular polynomials using MONomial Transformations:
//                      a C++ library for high precision integration of singular 
//                      polynomials of non-integer degree
//
// Authors:   Guido Lombardi, Davide Papapicco
//
// Institute: Politecnico di Torino
//            C.so Duca degli Abruzzi, 24 - Torino (TO), Italia
//            Department of Electronics and Telecommunications (DET)
//            Electromagnetic modelling and applications Research Group
//---------------------------------------------------------------------------------------

#include "Quasimont.h"

int main(int argc, char** argv)
{
  //  P_1(x) = ex^(e+1/4) + 5x^(pi/4) -x^(-1/2) +1 +10x^2
  std::vector<float50> coeff_sequence_bench_1 = {E, 5.0, -1.0, 1.0, 10.0};
  std::vector<float50> muntz_sequence_bench_1 = {E + 0.25, -PI/4, -0.5, 0, 2};

  quasimont(muntz_sequence_bench_1, coeff_sequence_bench_1);
  
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  std::cin.get();

  //  P_2(x) = x^(-e/3)
  std::vector<float50> coeff_sequence_bench_2 = {1};
  std::vector<float50> muntz_sequence_bench_2 = {-E/3};

  quasimont(muntz_sequence_bench_2, coeff_sequence_bench_2);
  
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  std::cin.get();

  //  P_3(x) = x^17 + x^35
  std::vector<float50> coeff_sequence_bench_3 = {1,1};
  std::vector<float50> muntz_sequence_bench_3 = {17,35};

  quasimont(muntz_sequence_bench_3, coeff_sequence_bench_3);
  
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  std::cin.get();

  return 0;
}