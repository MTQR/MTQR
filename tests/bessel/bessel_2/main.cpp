//---------------------------------------------------------------------------------------
// File:      tests/bessel/bessel_2/main.cpp
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

#include <gsl/gsl_sf_bessel.h>
#include "mtqr.h"

int main(int argc, char** argv){

  // Define the number of terms in the truncated series
  int N = 8;
  // Define the value of the exact definite integral obtained through symbolic computations in Wolfram Mathematica
  double I = 2.061676536411596;
  // Iterate over the number of terms in the truncated series to be added to the base case
  for(int k=0; k<=N; k++){
  	// Initialise the quadrature approximation to a null value
  	double In = 0.0;
  	// Compute the minimum and maximum monomial order of the Bessel function at iteration k
  	double nu_min = -0.5;
  	double nu_max = 0.333333333333333 + 2*k;
  	// Use those to apply the monomial transformation and obtain the new, optimised quadrature rule's nodes and weights
  	std::vector<std::vector<double>> quad_params = mtqr(nu_min, nu_max);
  	std::vector<double> nodes = quad_params.at(0);
  	std::vector<double> weights = quad_params.at(1);
  	// Compute the definite integral numerically using the MTQR with new nodes and weights
  	for(int j=0; j<nodes.size(); j++){
  		// Compute the value of the Bessel function of (fractional order) nu=nu_min when evaluated @ the j-th node of the loop
  		double J_nu_1 = gsl_sf_bessel_Jnu(nu_min, nodes[j]);
      double J_nu_2 = gsl_sf_bessel_Jnu(0.333333333333333, nodes[j]);
  		// Update the quadrature approximation with the current contribution
  		In += weights[j]*(J_nu_1 + J_nu_2);
  	}
  	// Print out the value of the quadrature at the current iteration of number of truncated terms
  	std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << "Iteration n° " << k << ":  In=" << In;
  	// Print out the value of the relative error
  	std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << ", En= " << fabs(I-In)/fabs(I);
  	// Print out the number of quadrature samples
  	std::cout << "   (n°samples=" << nodes.size() << ")" << "\n";
  }

  return 0;
}
