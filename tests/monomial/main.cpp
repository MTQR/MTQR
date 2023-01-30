//---------------------------------------------------------------------------------------
// File:      tests/monomial/main.cpp
//
// Library:   MTQR - Monomial Transformation Quadrature Rule:
//                   a C++ library for high-precision integration of 
//                   generalised polynomials of non-integer degree
//
// Authors:   Guido Lombardi, Davide Papapicco
//
// Institute: Politecnico di Torino
//            C.so Duca degli Abruzzi, 24 - Torino (TO), Italia
//            Department of Electronics and Telecommunications (DET)
//            Electromagnetic modelling and applications Research Group
//---------------------------------------------------------------------------------------

#include "mtqr.h"

int main(int argc, char** argv)
{
  //  P_1(x) = 0x^-0.5 + x^14 + 0x^28
  std::vector<float128> coeff_sequence = {0.0, 1.0, 0.0};
  std::vector<float128> muntz_sequence = {-0.5, 14, 28};

  mtqr(muntz_sequence, coeff_sequence);

  return 0;
}
