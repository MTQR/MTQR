//---------------------------------------------------------------------------------------
// File:      utilities/Utils.cpp
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
std::vector<type> castVector(const std::vector<float50>& input_vector, const type& type_infer)
{
  std::vector<type> output_vector;
  for(int k=0; k < input_vector.size(); k++)
  {
    output_vector.push_back(static_cast<type>(input_vector[k]));
  }
  return output_vector;
}
template std::vector<float50> castVector<float50>(const std::vector<float50>& input_vector, const float50& type_infer);
template std::vector<float128> castVector<float128>(const std::vector<float50>& input_vector, const float128& type_infer);
template std::vector<double> castVector<double>(const std::vector<float50>& input_vector, const double& type_infer);


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
float50 orderedInnerProduct(const std::vector<float50>& f_values, const std::vector<type>& weights)
{

  // Create vector storing the values of the projection of one input vector onto the other
  std::vector<float50> projection;
  // Compute the single terms of the projection
  for(int j=0; j < f_values.size(); j++)
  {
    projection.push_back(f_values[j]*static_cast<float50>(weights[j]));
  }
  // Sort the projection vector in ascending order
  sort(projection.begin(), projection.end());
  float50 inner_product = static_cast<float50>(0.0);

  // Compute the inner product by summing the ordered terms of the projection vector
  for(int k=0; k < projection.size(); k++)
  {
    inner_product += projection[k];
  }
  return inner_product;
}
template float50 orderedInnerProduct(const std::vector<float50>& f_values, const std::vector<float50>& weights);
template float50 orderedInnerProduct(const std::vector<float50>& f_values, const std::vector<float128>& weights);
template float50 orderedInnerProduct(const std::vector<float50>& f_values, const std::vector<double>& weights);


/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: equispaced_nodes = linspace(x_0, x_m, m)
//                
//          INPUT: - x_0 = starting node (infimum)
//                 - x_m = ending node (supremum)
//                 - m = number of sub-intervals
//
//         OUTPUT: - equispaced_nodes = array of m+1 equispaced <type> between x_0 and x_m 
//
//    DESCRIPTION: this method implements the linspace MATLAB function, generating a 
//                 linearly-spaced vector of nodes between a starting and ending point.
//
/////////////////////////////////////////////////////////////////////////////////////////

template<typename type>
std::vector<type> linspace(const type& start_type, const type& end_type, const int& num_steps)
{

  std::vector<type> linspaced;

  type start = static_cast<type>(start_type);
  type end = static_cast<type>(end_type);
  type steps = static_cast<type>(num_steps);

  if (steps == 0)
  { 
    return linspaced;
  }

  if (steps == 1) 
  {
    linspaced.push_back(start);
    return linspaced;
  }

  type delta = (end - start) / (steps - 1);

  for(int k = 0; k < steps - 1; k++)
  {
    linspaced.push_back(start + delta*k);
  }

  linspaced.push_back(end);
  return linspaced;
}
template std::vector<float50> linspace<float50>(const float50& start_type, const float50& end_type, const int& num_steps); // Template mock instantiation for non-inline function
template std::vector<float128> linspace<float128>(const float128& start_type, const float128& end_type, const int& num_steps); // Template mock instantiation for non-inline function
template std::vector<double> linspace<double>(const double& start_type, const double& end_type, const int& num_steps); // Template mock instantiation for non-inline function