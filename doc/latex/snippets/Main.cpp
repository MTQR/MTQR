//---------------------------------------------------------------------------------------
// File:      MyApp/Main.cpp
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
  //  P(x) = ex^(e+1/4) + 5x^(pi/4) -x^(-1/2) +1 +10x^2
  std::vector<float128> coeff_sequence = {E, 5.0, -1.0, 1.0, 10.0};
  std::vector<float128> muntz_sequence = {E + 0.25, -PI/4, -0.5, 0, 2};

  quasimont(muntz_sequence, coeff_sequence);

  return 0;
}