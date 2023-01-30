//---------------------------------------------------------------------------------------
// File:      src/mtqr.cpp
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

// Global variable controlling the primary module mode of execution
bool loud_mode = true;

// LOUD MODE
template<typename T>
void mtqr(std::vector<T>& muntz_sequence, std::vector<T>& coeff_sequence)
{
  // Print initial message and selects user's inputs
  auto input_data = manageData(muntz_sequence, coeff_sequence);

  // Extract beta_min, beta_max and n_min
  auto monomial_data = streamMonMapData(std::get<0>(input_data));

  // Compute order of the monomial transformation
  double transf_order = computeMapOrder(std::get<1>(input_data), std::get<1>(monomial_data));

  // Compute the new nodes and weights of the Monomial Transformation Quadrature Rule
  auto quad_data = computeQuadParams(transf_order, std::get<0>(monomial_data));

  // Cast the quadrature parameter in the most optimised f.p. format possible
  optimiseData(quad_data, muntz_sequence, coeff_sequence);
}
template void mtqr<float128>(std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence);
template void mtqr<double>(std::vector<double>& muntz_sequence, std::vector<double>& coeff_sequence);

// SILENT MODE
std::vector<std::vector<double>> mtqr(double lambda_min, double lambda_max)
{
  // Deactivate terminal's and files' output
  loud_mode = false;

  // Initialise input parameters of the Monomial transformation quadrature rule
  std::vector<double> muntz_sequence = {lambda_min, lambda_max};
  std::vector<double> coeff_sequence = {1.0, 1.0};

  // Print initial message and selects user's inputs
  auto input_data = manageData(muntz_sequence, coeff_sequence);

  // Extract beta_min, beta_max and n_min
  auto monomial_data = streamMonMapData(std::get<0>(input_data));

  // Compute order of the monomial transformation
  double transf_order = computeMapOrder(std::get<1>(input_data), std::get<1>(monomial_data));

  // Compute the new nodes and weights of the Monomial Transformation Quadrature Rule
  auto quad_data = computeQuadParams(transf_order, std::get<0>(monomial_data));

  // Cast the quadrature parameter in the most optimised f.p. format possible
  optimiseData(quad_data, muntz_sequence, coeff_sequence);

  // Generate double-precise new nodes and weights and export them in memory as output
  std::vector<double> nodes = castVector(std::get<0>(quad_data), std::numeric_limits<double>::epsilon());
  std::vector<double> weights = castVector(std::get<1>(quad_data), std::numeric_limits<double>::epsilon());
  return std::vector<std::vector<double>> {nodes, weights};
}