//---------------------------------------------------------------------------------------
// File:      src/Quasimont.cpp
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

/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: quasimont(muntz_sequence, coeff_sequence)
//                
//        INPUT: - muntz_sequence = sequence of real exponents of the polynomial
//               - coeff_sequence = sequence of real coefficients of the polynomial
//
//       OUTPUT: - no outputs
//
//    DESCRIPTION: access point of the focal module of the library where all the
//                 primary methods are instantiated according to the user's input
//
/////////////////////////////////////////////////////////////////////////////////////////

template<typename type>
void quasimont(std::vector<type>& muntz_sequence, std::vector<type>& coeff_sequence)
{
  // PRINT INITIAL MESSAGE AND SELECTS USER INPUTS
  auto input_data = manageData(muntz_sequence, coeff_sequence);

  // EXTRACT N_MIN, BETA_MIN AND BETA_MAX
  auto monomial_data = streamMonMapData(std::get<0>(input_data));

  // COMPUTE THE MONOMIAL TRANSFORMATION ORDER
  double transf_order = computeMapOrder(std::get<1>(input_data), std::get<1>(monomial_data));

  // COMPUTE AND EXPORT THE NEW G-L NODES & WEIGHTS
  auto quad_data = computeParamsGl(transf_order, std::get<0>(monomial_data));

  // CONVERTS AND EXPORTS NEW NODES AND WEIGHTS IN THE MOST OPTIMISED FLOATING-POINT FORMAT POSSIBLE
  optimiseData(quad_data, muntz_sequence, coeff_sequence);
}
template void quasimont(std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence);
template void quasimont(std::vector<double>& muntz_sequence, std::vector<double>& coeff_sequence);