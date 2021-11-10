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
//       FUNCTION: castVector(input_vector, output_vector)
//                
//          INPUT: - input_vector = vector of length n of type float50
//                 - input_vector = vector of length n of type T
//
//         OUTPUT: no outputs
//
//    DESCRIPTION: this method casts the float50 input vector to the an output vector
//                 with the same content but type T provided by the user.
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
//       FUNCTION: inner_product = ordereInnerProduct(input_vector)
//                
//          INPUT: - input_vector = input vector of type float1k
//
//         OUTPUT: - inner_product = 
//
//    DESCRIPTION: in multiple instances throughout the library, addition of terms close
//                 to the format epsilon is performed (especially when computing the
//                 quadrature). To avoid numeric cancellation of these terms we must 
//                 therefore sum those smallest values first in the highest precision
//                 possible. This method takes as input the terms of the inner product
//                 of the terms of two vectors of length n, sorts it in ascending order
//                 and sums along n, thereby assuring that no numerical cancellation
//                 occurs.
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