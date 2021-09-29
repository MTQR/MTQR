//------------------------------------------------------------------------------------------------------------------
// File:      utilities/Linspace.cpp
//
// Library:   QUASIMONT - QUAdrature of SIngular polynomials using MONomial Transformations:
//                        a C++ library for high precision integration of singular polynomials of non-integer degree
//
// Authors:   Davide Papapicco, Guido Lombardi, PhD
//
// Institute: Politecnico di Torino
//            C.so Duca degli Abruzzi, 24 - Torino (TO), Italia
//            Department of Electronics and Telecommunications (DET)
//            Electromagnetic modelling and applications Research Group
//------------------------------------------------------------------------------------------------------------------

#include "Quasimont.h"


//******************************************************************************************************************
//******************************************************************************************************************
//――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
//
//       FUNCTION: equispaced_nodes = linspace(x_0, x_m, m)
//                
//          INPUT: - x_0 = starting node (infimum)
//                 - x_m = ending node (supremum)
//                 - m = number of sub-intervals
//
//         OUTPUT: - equispaced_nodes = array of m+1 equispaced <type> between x_0 and x_m 
//
//    DESCRIPTION: this method implements the linspace MATLAB function, generating a linearly-spaced vector of nodes
//                 between a starting and ending point.


template<typename type>
void castVector(const std::vector<float50>& input_vector, std::vector<type>& output_vector)
{

  for(int k=0; k < input_vector.size(); k++)
  {
    output_vector.push_back(static_cast<type>(input_vector[k]));
  }
}
template void castVector<float50>(const std::vector<float50>& input_vector, std::vector<float50>& output_vector);
template void castVector<float128>(const std::vector<float50>& input_vector, std::vector<float128>& output_vector);
template void castVector<double>(const std::vector<float50>& input_vector, std::vector<double>& output_vector);
//――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――


//******************************************************************************************************************
//******************************************************************************************************************
//――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
//
//       FUNCTION: equispaced_nodes = linspace(x_0, x_m, m)
//                
//          INPUT: - x_0 = starting node (infimum)
//                 - x_m = ending node (supremum)
//                 - m = number of sub-intervals
//
//         OUTPUT: - equispaced_nodes = array of m+1 equispaced <type> between x_0 and x_m 
//
//    DESCRIPTION: this method implements the linspace MATLAB function, generating a linearly-spaced vector of nodes
//                 between a starting and ending point.


float1k orderedInnerProduct(std::vector<float1k>& input_vector)
{

  float1k inner_product = 0;
  sort(input_vector.begin(), input_vector.end());

  for(int k=0; k < input_vector.size(); k++)
  {

    inner_product += input_vector[k];
  }

  return inner_product;
}


//******************************************************************************************************************
//******************************************************************************************************************
//――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
//
//       FUNCTION: equispaced_nodes = linspace(x_0, x_m, m)
//                
//          INPUT: - x_0 = starting node (infimum)
//			           - x_m = ending node (supremum)
//                 - m = number of sub-intervals
//
//		     OUTPUT: - equispaced_nodes = array of m+1 equispaced <type> between x_0 and x_m 
//
//    DESCRIPTION: this method implements the linspace MATLAB function, generating a linearly-spaced vector of nodes
//                 between a starting and ending point.


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
template std::vector<float128> linspace<float128>(const float128& start_type, const float128& end_type, const int& num_steps); // Template mock instantiation for non-inline function
//――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――