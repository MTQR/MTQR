//---------------------------------------------------------------------------------------
// File:      src/Quasimont.cpp
//
// Library:   QUASIMONT-QUAdrature of SIngular polynomials using MONomial Transformations:
//                      a C++ library for high precision integration of generalised 
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

// Global variable controlling QUASIMONT's terminal behaviour
bool loud_mode = true;

/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: quasimont(muntz_sequence, coeff_sequence)
//                
//          INPUT: - muntz_sequence = sequence of real exponents of the polynomial
//                 - coeff_sequence = sequence of real coefficients of the polynomial
//
//         OUTPUT: - no outputs
//
//    DESCRIPTION: access point of the primary module of the library where all the
//                 methods concerning the computation of the monomial transformation 
//                 quadrature rule's params are instantiated according to user's input.
//
/////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
void quasimont(std::vector<T>& muntz_sequence, std::vector<T>& coeff_sequence)
{
  // PRINT INITIAL MESSAGE AND SELECTS USER INPUTS
  auto input_data = manageData(muntz_sequence, coeff_sequence);

  // EXTRACT N_MIN, BETA_MIN AND BETA_MAX
  auto monomial_data = streamMonMapData(std::get<0>(input_data));

  // COMPUTE THE MONOMIAL TRANSFORMATION ORDER
  double transf_order = computeMapOrder(std::get<1>(input_data), std::get<1>(monomial_data));

  // COMPUTE AND EXPORT THE NEW NODES & WEIGHTS OF THE MONOMIAL TRANSFORMATION QUADRATURE RULE
  auto quad_data = computeQuadParams(transf_order, std::get<0>(monomial_data));

  // CONVERTS AND EXPORTS NEW NODES AND WEIGHTS IN THE MOST OPTIMISED FLOATING-POINT FORMAT POSSIBLE
  optimiseData(quad_data, muntz_sequence, coeff_sequence);
}
template void quasimont<float128>(std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence);
template void quasimont<double>(std::vector<double>& muntz_sequence, std::vector<double>& coeff_sequence);

// Silent mode overloading
std::vector<std::vector<double>> quasimont(double lambda_min, double lambda_max)
{
  // ACTIVATE TERMINAL'S OUTPUT SUPPRESSION FLAG
  loud_mode = false;

  // MAKE INPUT PARAMETERS FOR QUASIMONT'S WORK-FLOW
  std::vector<double> muntz_sequence = {lambda_min, lambda_max};
  std::vector<double> coeff_sequence = {1.0, 1.0};

  // PRINT INITIAL MESSAGE AND SELECTS USER INPUTS
  auto input_data = manageData(muntz_sequence, coeff_sequence);

  // EXTRACT N_MIN, BETA_MIN AND BETA_MAX
  auto monomial_data = streamMonMapData(std::get<0>(input_data));

  // COMPUTE THE MONOMIAL TRANSFORMATION ORDER
  double transf_order = computeMapOrder(std::get<1>(input_data), std::get<1>(monomial_data));

  // COMPUTE AND EXPORT THE NEW NODES & WEIGHTS OF THE MONOMIAL TRANSFORMATION QUADRATURE RULE
  auto quad_data = computeQuadParams(transf_order, std::get<0>(monomial_data));

  // EXPORTS NEW NODES AND WEIGHTS IN THE MOST OPTIMISED FLOATING-POINT FORMAT POSSIBLE
  optimiseData(quad_data, muntz_sequence, coeff_sequence);

  // CONVERTS THE NEW NODES AND WEIGHTS IN DOUBLE F.P. FORMAT
  std::vector<double> nodes = castVector(std::get<0>(quad_data), std::numeric_limits<double>::epsilon());
  std::vector<double> weights = castVector(std::get<1>(quad_data), std::numeric_limits<double>::epsilon());
  
  return std::vector<std::vector<double>> {nodes, weights};
}