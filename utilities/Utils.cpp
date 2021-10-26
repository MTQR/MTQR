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

template<typename inputType, typename outputType>
void castVector(const std::vector<inputType>& input_vector, std::vector<outputType>& output_vector)
{
  for(int k=0; k < input_vector.size(); k++)
  {
    output_vector.push_back(static_cast<outputType>(input_vector[k]));
  }
}
template void castVector<float50, float50>(const std::vector<float50>& input_vector, std::vector<float50>& output_vector);
template void castVector<float50, float128>(const std::vector<float50>& input_vector, std::vector<float128>& output_vector);
template void castVector<float50, double>(const std::vector<float50>& input_vector, std::vector<double>& output_vector);
template void castVector<float128, float50>(const std::vector<float128>& input_vector, std::vector<float50>& output_vector);
template void castVector<float128, float128>(const std::vector<float128>& input_vector, std::vector<float128>& output_vector);
template void castVector<float128, double>(const std::vector<float128>& input_vector, std::vector<double>& output_vector);
template void castVector<double, float50>(const std::vector<double>& input_vector, std::vector<float50>& output_vector);
template void castVector<double, float128>(const std::vector<double>& input_vector, std::vector<float128>& output_vector);
template void castVector<double, double>(const std::vector<double>& input_vector, std::vector<double>& output_vector);


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
type orderedInnerProduct(const std::vector<type>& f_values, const std::vector<type>& weights)
{
  //std::cout << "\n ** Element-wise Product";
  // Create vector storing the values of the projection of one input vector onto the other
  std::vector<type> projection;
  // Compute the single terms of the projection
  for(int j=0; j < f_values.size(); j++)
  {
    projection.push_back(f_values[j]*weights[j]);
    /*std::cout << "\n **** f(x_" << j << ") = " << f_values[j];
    std::cout << "\n **** w_" << j << ") = " << weights[j];
    std::cout << "\n **** f(x_" << j << ")*w_" << j << " = " << f_values[j]*weights[j];*/
  }
  // Sort the projection vector in ascending order
  sort(projection.begin(), projection.end());

  //std::cout << "\n ** Ordered Inner Product";
  type inner_product = 0;
  // Compute the inner product by summing the ordered terms of the projection vector
  for(int k=0; k < projection.size(); k++)
  {
    inner_product += projection[k];
    //std::cout << "\n **** sum = " << inner_product;
  }
  return inner_product;
}
template float50 orderedInnerProduct(const std::vector<float50>& f_values, const std::vector<float50>& weights);
template float128 orderedInnerProduct(const std::vector<float128>& f_values, const std::vector<float128>& weights);
template double orderedInnerProduct(const std::vector<double>& f_values, const std::vector<double>& weights);


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






void plot(const int& num_nodes)//, const double& min_exponent, const double& max_exponent)
{

  //int print_startpoint = static_cast<int>(ceil(min_exponent)) -10;
  int print_startpoint = 0;
  int print_endopoint = 200;// static_cast<int>(ceil(max_exponent)) + 10;
  double dp_epsilon = std::numeric_limits<double>::epsilon();

  std::ofstream epsilon;

  std::string epsilon_file = "./output/Epsilon.csv";
  epsilon.open(epsilon_file.c_str());

  epsilon << std::setprecision(std::numeric_limits<float1k>::max_digits10) << print_startpoint << "," << dp_epsilon << "\n" << print_endopoint << "," << dp_epsilon;
  epsilon.close();

  std::ofstream f_exact, f_exact_env;

  f_exact.open("./output/ExactError.csv");
  f_exact_env.open("./output/ExactErrorEnv.csv");

  int n_samples = 13517;
  std::vector<float128> lambda = linspace(static_cast<float128>(-1), static_cast<float128>(200), n_samples);

  std::string envelope = "yes";
  std::string no_envelope = "no";

  std::cout << "\nComputing the error estimates with " << num_nodes << " quadrature nodes:  " << std::endl;

  for(int k = 0; k < n_samples - 2; k++)
  {
    
    f_exact << std::setprecision(std::numeric_limits<float1k>::max_digits10) << lambda[k+1] << "," << computeErrorEstimate(lambda[k+1], num_nodes, no_envelope) << "\n";
    f_exact_env << std::setprecision(std::numeric_limits<float1k>::max_digits10) << lambda[k+1] << "," << computeErrorEstimate(lambda[k+1], num_nodes, envelope) << "\n";

    printProgressBar(k, n_samples - 2);

  }

  epsilon.close();
  f_exact.close();
  f_exact_env.close();

  std::cout << "\nComputation completed!" << std::endl;
}

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