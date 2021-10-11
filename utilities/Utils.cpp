//---------------------------------------------------------------------------------------
// File:      utilities/Utils.cpp
//
// Library:   QUASIMONT-QUAdrature of SIngular polynomials using MONomial Transformations:
//                      a C++ library for high precision integration of singular 
//                      polynomials of non-integer degree
//
// Authors:   Guido Lombardi, PhD, Davide Papapicco
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
template std::vector<float128> linspace<float128>(const float128& start_type, const float128& end_type, const int& num_steps); // Template mock instantiation for non-inline function

/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: printProgressBar(iterator, number_of_iterations)
//                
//          INPUT: - iterator = value of the incremending integer index in the 'for' loop
//                 - number_of_iterations = end-value of the iterator in the 'for' loop
//          
//         OUTPUT: no outputs
//
//    DESCRIPTION: in a traditional for(int k=0; k<num_k; k++) this function prints a 
//                 progress bar on the terminal line to represent the status of
//                 completion of the loop.
//
/////////////////////////////////////////////////////////////////////////////////////////

void printProgressBar(const int& iter, const int& num_iter)
{

  int bar_length = 70;

  float ratio = static_cast<float>(iter)/static_cast<float>(num_iter);
  int progress = ceil(ratio*70);
  int progress_prctg = ceil(ratio*100.0);

  std::cout << "[";
  
  for(int k = 0; k <= bar_length; k++)
  {
    if(k <= progress)
    {
      std::cout << "#";
    }
    else
    {
      std::cout << ".";
    }
  }

  std::cout << "] " << std::setprecision(2) << progress_prctg << " %\r";
  std::cout.flush();
}