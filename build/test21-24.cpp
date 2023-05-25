#include <gsl/gsl_sf_bessel.h>
#include <chrono>
#include "mtqr.h"

// clear && rm outputs.txt && g++ -w example1-4.cpp -o example1-4 -I../include -L. -lmtqr -lm -lquadmath -lgsl -lgslcblas && ./example1-4

int main(int argc, char** argv){
  
  std::vector<double> In;
  std::vector<double> Tn;

  std::vector<std::vector<double>> nu_s = {{0.0, 1.0}, {-1.0/3.0, 0.0}, {-1.0/2.0, -1.0/3.0}, {0.0, -3.141592653589793/4.0}};
  std::vector<std::vector<double>> coef = {{1.0, 3.0/2.0}, {1.0, 3.0}, {1, 1}, {1, 1}};
  std::vector<double> m = {-1.0/2.0, 1.0/6.0, 0, 0};

  double integrand;
  std::vector<std::vector<double>> lambdas = {{1.0/2.0, 17.0/2.0}, {-1.0/6.0, 71.0/6.0}, {-5.0/6.0, 43.0/6.0}, {-3.141592653589793/4.0, 8.0 - 3.141592653589793/4.0}};
  
  std::vector<std::vector<double>> quad_params;
  std::vector<double> nodes, weights;

  for(int k=0; k<m.size(); k++){
    // Initialise dummy variable to hold the value of the approximated integral
    double dummy = 0.0;
  	// Get the two coefficients of the argument of the integrand products of Bessel functions at iteration k
  	double coeff_1 = coef[k][0];
  	double coeff_2 = coef[k][1];
  	// Get the two orders of the integrand products of Bessel functions at iteration k
  	double nu_1 = nu_s[k][0];
  	double nu_2 = nu_s[k][1];
    // Get the minimum and maximum monomial order of the polynomial series approximating the integrands in [0,1]
    double lambda_min = lambdas[k][0];
    double lambda_max = lambdas[k][1];
  	
    // Use those to apply the monomial transformation and obtain the new, optimised quadrature rule's nodes and weights
  	auto start = std::chrono::high_resolution_clock::now();
    quad_params = mtqr(lambda_min, lambda_max);
    auto stop = std::chrono::high_resolution_clock::now();

  	nodes = quad_params.at(0);
  	weights = quad_params.at(1);

  	// Debug
  	// std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << "k=" << k+1 << ", nu_1=" << nu_1 << ", nu_2=" << nu_2 << ", a_1=" << coeff_1 << ", a_2=" << coeff_2 << "\n";
  	
    // Compute the definite integral numerically using the MTQR with new nodes and weights
  	for(int j=0; j<nodes.size(); j++){
  		// Compute the value of the Bessel function of (fractional order) nu=nu_min when evaluated @ the j-th node of the loop
  		integrand = pow(nodes[j], m[k])*gsl_sf_bessel_Jnu(nu_1, coeff_1*nodes[j])*gsl_sf_bessel_Jnu(nu_2, coeff_2*nodes[j]);
  		// Update the quadrature approximation with the current contribution
  		dummy += weights[j]*integrand;
  		// Debug
  		/*
  		std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
  				  << "  j=" << j+1
  				  << ", x=" << nodes[j]
  				  << ", w=" << weights[j]
  				  << "  ==>  " << integrand
  				  << ", In=" << dummy
  				  << "\n";
  		*/
  	}

    // Update the array of quadratures with the computed value held in the dummy variable
    In.push_back(dummy);

    // Compute and update the computational time
    auto time = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    Tn.push_back(time.count()*1e-6);
  }

  // Print the results
  for(int k=0; k<In.size(); k++){
    std::cout << "k=" << k+1 << ",   In = " << std::setprecision(std::numeric_limits<double>::max_digits10) << In[k] 
    					     << ",   Tn = " << std::setprecision(std::numeric_limits<float>::max_digits10) << Tn[k] << " seconds\n";
  }

  return 0;
}