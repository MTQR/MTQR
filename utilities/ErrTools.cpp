//---------------------------------------------------------------------------------------
// File:      utilities/ErrTools.cpp
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
#include <boost/multiprecision/mpfr.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/sin_pi.hpp>

typedef boomp::number<boomp::mpfr_float_backend<50>> float50; // Definition of f.p. datatype with 50 decimal digits of precision

using boost::math::policies::domain_error;
using boost::math::policies::errno_on_error;

typedef boost::math::policies::policy<
   boost::math::policies::domain_error<boost::math::policies::errno_on_error>,
   boost::math::policies::domain_error<boost::math::policies::errno_on_error>,
   boost::math::policies::domain_error<boost::math::policies::errno_on_error>,
   boost::math::policies::domain_error<boost::math::policies::errno_on_error>
> gamma_policy; // Sets the domain error handling policy for the tgamma boost special function

/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: equispaced_nodes = linspacedVector(x_0, x_m, m)
//                
//          INPUT: - x_0 = starting node (infimum)
//                 - x_m = ending node (supremum)
//                 - m = number of sub-intervals
//
//         OUTPUT: - equispaced_nodes = array of m+1 equispaced float50 nodes between
//                                      x_0 and x_m 
//
//    DESCRIPTION: this method implements MATLAB's 'linspace' function, generating a 
//                 linearly-spaced vector of nodes between a starting and ending point.
//
/////////////////////////////////////////////////////////////////////////////////////////

std::vector<float50> linspacedVector(const float50& start_type, const float50& end_type, const int& num_steps)
{
  std::vector<float50> linspaced;

  float50 start = static_cast<float50>(start_type);
  float50 end = static_cast<float50>(end_type);
  float50 steps = static_cast<float50>(num_steps);

  if (steps == 0)
  { 
    return linspaced;
  }

  if (steps == 1) 
  {
    linspaced.push_back(start);
    return linspaced;
  }

  float50 delta = (end - start) / (steps - 1);
  for(int k = 0; k < steps - 1; k++)
  {
    linspaced.push_back(start + delta*k);
  }
  linspaced.push_back(end);
  return linspaced;
}

/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: En = aPrioriAsympEstimate(lambda, n)
//                
//          INPUT: - lambda = real value of the non-integer degree of the monomial
//                 - n = number of quadrature nodes
//
//         OUTPUT: - En = exact asymptotic estimate of the G-L quadrature error 
//                        computed 1.9 and 1.8
//
//    DESCRIPTION: as a demonstration, this routine computes the values of the error
//                 estimate for a given value of n and lambda.
//
/////////////////////////////////////////////////////////////////////////////////////////

template<typename type>
type aPrioriAsympEstimate(const type& lambda, const int& num_nodes)
{
  type common_factor = 2*(pow(static_cast<type>(2.0),-(static_cast<type>(1.0)+lambda)));
  type denominator = 2*num_nodes + lambda;

  type b1 = (boost::math::tgamma(static_cast<type>(2.0)*lambda, gamma_policy())*boost::math::tgamma(static_cast<type>(2.0*num_nodes) - lambda, gamma_policy()))/boost::math::tgamma(static_cast<type>(2.0)*lambda + static_cast<type>(2.0*num_nodes) - lambda, gamma_policy());
  type b2 = (boost::math::tgamma(static_cast<type>(2.0)*lambda, gamma_policy())*boost::math::tgamma(static_cast<type>(2.0) + static_cast<type>(2.0*num_nodes) - lambda, gamma_policy()))/boost::math::tgamma(static_cast<type>(2.0)*lambda + static_cast<type>(2.0) + static_cast<type>(2.0*num_nodes) - lambda, gamma_policy());

  type exact = common_factor*(pow(static_cast<type>(2.0), -lambda))*lambda*((b1/denominator) - b2/(2 + denominator));
  
  if(lambda + static_cast<type>(0.5) - 2*num_nodes >= 0)
  {
    return fabs(exact*boost::math::sin_pi(lambda));
  }
  else
  {
    return fabs(exact);
  }
}
template float50 aPrioriAsympEstimate<float50>(const float50& input_lambda, const int& num_nodes);

/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: plot(num_nodes, beta_min, beta_max)
//                
//          INPUT: - num_nodes = number of quadrature nodes
//                 - beta_min = minimum value of the exponents 
//                 - beta_max = maximum value of the exponents
//          
//         OUTPUT: no outputs
//
//    DESCRIPTION: this method plots the asymptotic behaviour of the error estimate 
//                 computed in function 'aPrioriAsympEstimate' and for a fixed value of n.
//                 Several examples of its plots are superimposed in Figure 1.1 of the
//                 doc/UserManual.pdf. 
//
/////////////////////////////////////////////////////////////////////////////////////////

template<typename type>
void plot(const int& num_nodes, const type& beta_min, const type& beta_max)
{
  std::cout << "\n\n *** Plotting the a-priori asymptotic error estimate with N=" << num_nodes << " ***" << std::endl;
  
  int print_startpoint = -1;
  int print_endopoint = static_cast<int>(ceil(beta_max)) + 10;

  std::ofstream epsilon;
  std::string epsilon_file = "./estimate/Epsilon.csv";
  epsilon.open(epsilon_file.c_str());
  epsilon << std::setprecision(std::numeric_limits<type>::max_digits10) << print_startpoint << "," << EPS << "\n" << print_endopoint << "," << EPS;
  epsilon.close();

  std::string csv = ".csv";
  std::string error_file = "./estimate/EstimateN=" + std::to_string(num_nodes) + csv;
  std::ofstream error;
  error.open(error_file);

  int n_samples = (beta_max - beta_min)*100;
  std::vector<type> lambda = linspacedVector(static_cast<type>(print_startpoint), static_cast<type>(print_endopoint), n_samples);

  for(int k = 0; k < n_samples - 2; k++)
  {
    error << std::setprecision(std::numeric_limits<type>::max_digits10) << lambda[k+1] << "," << aPrioriAsympEstimate(lambda[k+1], num_nodes) << "\n";
    printProgressBar(k, n_samples - 2);
  }

  std::string plot_file = "bash ./plot.sh ";
  std::string cmd = plot_file + std::to_string(num_nodes) + " " + std::to_string(print_startpoint) + " " + std::to_string(print_endopoint);
  system(cmd.c_str());

  error.close();
}
template void plot<float50>(const int& num_nodes, const float50& beta_min, const float50& beta_max);

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

/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: main(argc, argv)
//
//    DESCRIPTION: a standard C++ main function dictating the sequence of methods to be
//                 executed by the secondary module of the library.
//
/////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  int plot_nodes = 42;

  system("mkdir -p estimate");
  std::ofstream data("./estimate/TabBetas.csv");

  int num_samples = 15000, prev_j=0;
  std::vector<float50> lambda = linspacedVector(static_cast<float50>(-1.0), static_cast<float50>(1300.0), num_samples);

  float50 exact, tab_beta_min, tab_beta_max;
  int min_nodes = 5, max_nodes = 50;

  for(int n = min_nodes; n <= max_nodes; n++)
  {
    int nodes = 2*n;
    std::cout << "\nComputing for n = " << nodes << std::endl;

    for(int k = 0; k < num_samples - 2; k++)
    {
      std::cout << "    beta_min = " << lambda[k+1] << "\r";

      exact = aPrioriAsympEstimate(lambda[k+1], nodes);
      if(exact <= EPS)
      {
        std::cout << std::endl;
        tab_beta_min = lambda[k+1];

        if(prev_j < k+1)
        {
          prev_j = k+1;
        }
        for(int j = prev_j+1; j < num_samples - 2; j++)
        {
          std::cout << "    beta_max = " << lambda[j] << "\r";

          exact = aPrioriAsympEstimate(lambda[j], nodes);
          if(exact >= EPS)
          {
            tab_beta_max = lambda[j];
            prev_j = j;

            break;
          }
        }
        break;
      }
    }
    data << std::setprecision(std::numeric_limits<float50>::max_digits10) << nodes << "," << tab_beta_min << "," << tab_beta_max << "\n";

    if(nodes==plot_nodes)
    {
      plot(nodes, tab_beta_min, tab_beta_max);
      return 0;
    }
  }
  data.close();

  std::cout << "\n\n ** Data available in the 'estimate' subdirectory **\n ** Done! **\n" << std::endl;;
}