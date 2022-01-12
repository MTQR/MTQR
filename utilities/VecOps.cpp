//---------------------------------------------------------------------------------------
// File:      utilities/VecOps.cpp
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
//       FUNCTION: output_vector = castVector(input_vector, type_output)
//                
//          INPUT: - input_vector = vector of length n of type of type T1
//                 - type_output = floating-point data-type T2
//
//         OUTPUT: - output_vector = input_vector casted in type T2
//
//    DESCRIPTION: this method casts the T1=float128 input vector to the an output vector
//                 with the same content but type T2 specified by the user.
//
/////////////////////////////////////////////////////////////////////////////////////////

template<typename type>
std::vector<type> castVector(const std::vector<float128>& input_vector, const type& type_infer)
{
  std::vector<type> output_vector;
  for(int k=0; k < input_vector.size(); k++)
  {
    output_vector.push_back(static_cast<type>(input_vector[k]));
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
//                 arranged in ascending order and sums them thereby assuring that no 
//                 numerical cancellation occurs.
//
/////////////////////////////////////////////////////////////////////////////////////////

template<typename type>
float128 doubleDotProduct(const std::vector<float128>& f_values, const std::vector<type>& weights)
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