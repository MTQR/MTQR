#include <gsl/gsl_sf_bessel.h>
#include <chrono>
#include "mtqr.h"

// clear && rm outputs.txt && g++ -w comparison.cpp -o comparison -I../include -L. -lmtqr -lm -lquadmath -lgsl -lgslcblas && ./comparison

int main(int argc, char** argv){
  
  std::vector<double> In;
  std::vector<double> En;
  std::vector<double> Tn;

  std::vector<double> I = {0.20723625024316799, 0.00011234154868471322, 1.4737950687864697, 1.3863773927287198, 4.1966646744359134};
  std::vector<std::vector<double>> nu_s = {{0.0, 1.0}, {0.0, 5.0}, {0.0, 1.0/2.0}, {0.0, -1.0/2.0}, {-1.0/2.0, -1.0/3.0}};
  std::vector<int> m = {0, -4, -1, 0, 0};

  double integrand;
  std::vector<std::vector<double>> lambdas = {{1.0, 9.0}, {1.0, 9.0}, {-1.0/2.0, 15.0/2.0}, {-1.0/2.0, 15.0/2.0}, {-5.0/6.0, 43.0/6.0}};
  
  std::vector<std::vector<double>> quad_params;
  std::vector<double> nodes, weights;

  for(int k=0; k<I.size(); k++){
    // Initialise dummy variable to hold the value of the approximated integral
    double dummy = 0.0;
  	// Get the two orders of the integrand products of Bessel functions at iteration k
  	double nu_1 = nu_s[k][1];
  	double nu_2 = nu_s[k][2];

    // get the minimum and maximum monomial order of the polynomial series approximating the integrands in [0,1]
    double lambda_min = lambdas[k][1];
    double lambda_max = lambdas[k][2];
  	
    // Use those to apply the monomial transformation and obtain the new, optimised quadrature rule's nodes and weights
  	auto start = std::chrono::high_resolution_clock::now();
    quad_params = mtqr(lambda_min, lambda_max);
    auto stop = std::chrono::high_resolution_clock::now();

  	nodes = quad_params.at(0);
  	weights = quad_params.at(1);
  	
    // Compute the definite integral numerically using the MTQR with new nodes and weights
  	for(int j=0; j<nodes.size(); j++){
  		// Compute the value of the Bessel function of (fractional order) nu=nu_min when evaluated @ the j-th node of the loop
  		integrand = pow(nodes[j], m[k])*gsl_sf_bessel_Jnu(nu_1, nodes[j])*gsl_sf_bessel_Jnu(nu_2, nodes[j]);
  		// Update the quadrature approximation with the current contribution
  		dummy += weights[j]*integrand;
  	}

    // Update the array of quadratures with the computed value held in the dummy variable
    In.push_back(dummy);

    // Compute and update the relative error
    En.push_back(abs(In[k]-I[k])/abs(I[k]));

    // Compute and update the computational time
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    Tn.push_back(time.count());
  }

  // Print the results
  for(int k=0; k<In.size(); k++){
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << "k=" << k << ",   En = " << En[k] << ",   Tn = " << Tn[k] << "*10-3 seconds\n";
  }

  return 0;
}