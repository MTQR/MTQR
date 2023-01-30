//---------------------------------------------------------------------------------------
// File:      utilities/vector_ops.cpp
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

/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: output_vector = castVector(input_vector, type_output)
//                
//          INPUT: - input_vector = vector of length n of type of type float128
//                 - type_output = floating-point templatised type T
//
//         OUTPUT: - output_vector = input_vector casted in type T
//
//    DESCRIPTION: this method casts the float128 input vector to the an output vector
//                 with the same content but type T specified by the user.
//
/////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
std::vector<T> castVector(const std::vector<float128>& input_vector, const T& type_infer)
{
  std::vector<T> output_vector;
  for(int k=0; k < input_vector.size(); k++)
  {
    output_vector.push_back(static_cast<T>(input_vector[k]));
  }
  return output_vector;
}
template std::vector<float128> castVector<float128>(const std::vector<float128>& input_vector, const float128& type_infer);
template std::vector<double> castVector<double>(const std::vector<float128>& input_vector, const double& type_infer);

/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: inner_product = doubleDotproduct(x, w)
//                
//          INPUT: - x = vector containing the quadrature nodes in format float128
//                 - w = vector containing the quadrature weights in either float128 or 
//                       double
//
//         OUTPUT: - inner_product = precise inner product between x and w avoiding
//                                   numerical cancellation
//
//    DESCRIPTION: This method takes two input vectors of the same length n, sorts
//                 their element-wise multiplication in a new vector with cells
//                 arranged in ascending order and sums them thereby limiting numerical
//                 cancellation errors.
//
/////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
float128 doubleDotProduct(const std::vector<float128>& f_values, const std::vector<T>& weights)
{

  // Create vector storing the values of the projection of one input vector onto the other
  std::vector<float128> projection;
  // Compute the single terms of the projection
  for(int j=0; j < f_values.size(); j++)
  {
    projection.push_back(f_values[j]*static_cast<float128>(weights[j]));
  }
  // Sort the projection vector in ascending order
  sort(projection.begin(), projection.end());
  float128 inner_product = static_cast<float128>(0.0);

  // Compute the inner product by summing the ordered terms of the projection vector
  for(int k=0; k < projection.size(); k++)
  {
    inner_product += projection[k];
  }
  return inner_product;
}
template float128 doubleDotProduct(const std::vector<float128>& f_values, const std::vector<float128>& weights);
template float128 doubleDotProduct(const std::vector<float128>& f_values, const std::vector<double>& weights);