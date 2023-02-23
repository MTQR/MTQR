//---------------------------------------------------------------------------------------
// File:      tests/bessels/main.cpp
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
  int N=4;
  double I, In;
  double nu_min, nu_max;
  double J_nu, J_nu_1, J_nu_2;
  std::vector<double> nodes, weights;
  std::vector<std::vector<double>> quad_params;

  /*  BESSEL TEST 1  */
  // Define the value of the exact definite integral obtained through symbolic computations in Wolfram Mathematica
  I = 1.443411848585210;
  // Iterate over the number of terms in the truncated series to be added to the base case
  for(int k=1; k<=N; k++){
  	// Initialise the quadrature approximation to a null value
  	In = 0.0;
  	// Compute the minimum and maximum monomial order of the Bessel function at iteration k
  	nu_min = -0.5;
  	nu_max = -0.5 + 2*k;
  	// Use those to apply the monomial transformation and obtain the new, optimised quadrature rule's nodes and weights
  	quad_params = mtqr(nu_min, nu_max);
  	nodes = quad_params.at(0);
  	weights = quad_params.at(1);
  	// Compute the definite integral numerically using the MTQR with new nodes and weights
  	for(int j=0; j<nodes.size(); j++){
  		// Compute the value of the Bessel function of (fractional order) nu=nu_min when evaluated @ the j-th node of the loop
  		J_nu = gsl_sf_bessel_Jnu(nu_min, nodes[j]);
  		// Update the quadrature approximation with the current contribution
  		In += weights[j]*J_nu;
  	}
  	// Print out the value of the quadrature at the current iteration of number of truncated terms
  	std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << "Iteration n° " << k << ":  In=" << In;
  	// Print out the value of the relative error
  	std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << ", En=" << fabs(I-In)/fabs(I);
  	// Print out the number of quadrature samples
  	std::cout << "   (n°samples=" << nodes.size() << ")" << "\n";
  }
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  std::cin.get();

  /*  BESSEL TEST 2  */
  I = 2.061676536411596;
  for(int k=1; k<=N; k++){
    In = 0.0;
    nu_min = -0.5;
    nu_max = 0.333333333333333 + 2*k;
    
    quad_params = mtqr(nu_min, nu_max);
    nodes = quad_params.at(0);
    weights = quad_params.at(1);
    for(int j=0; j<nodes.size(); j++){
      J_nu_1 = gsl_sf_bessel_Jnu(nu_min, nodes[j]);
      J_nu_2 = gsl_sf_bessel_Jnu(0.333333333333333, nodes[j]);

      In += weights[j]*(J_nu_1 + J_nu_2);
    }
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << "Iteration n° " << k << ":  In=" << In;
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << ", En= " << fabs(I-In)/fabs(I);
    std::cout << "   (n°samples=" << nodes.size() << ")" << "\n";
  }
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  std::cin.get();

  /*  BESSEL TEST 3  */
  I = 0.2021793666591346;
  for(int k=1; k<=N; k++){In = 0.0;
    nu_min = 0.5;
    nu_max = 0.5 + 4*k;
    
    quad_params = mtqr(nu_min, nu_max);
    nodes = quad_params.at(0);
    weights = quad_params.at(1);
    for(int j=0; j<nodes.size(); j++){
      J_nu_1 = gsl_sf_bessel_Jnu(-nu_min, nodes[j]);
      J_nu_2 = gsl_sf_bessel_Jnu(1.0, nodes[j]);

      In += weights[j]*J_nu_1*J_nu_2;
    }
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << "Iteration n° " << k << ":  In=" << In;
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << ", En= " << fabs(I-In)/fabs(I);
    std::cout << "   (n°samples=" << nodes.size() << ")" << "\n";
  }
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  std::cin.get();

  return 0;
}
