#include <gsl/gsl_sf_bessel.h>
#include <chrono>
#include "mtqr.h"

// clear && rm outputs.txt && g++ -w new_test25.cpp -o test25 -I../include -L. -lmtqr -lm -lquadmath -lgsl -lgslcblas && ./test25

int main(int argc, char** argv){
  
  std::vector<double> In;
  std::vector<double> Tn;

  std::vector<double> mu = {-0.8456343634308191, -0.7516749897162837, -0.6577156160017482, -0.5637562422872127, -0.4697968685726773, -0.3758374948581418, -0.2818781211436064, -0.1879187474290709, -0.09395937371453546, 0.0};
  std::vector<double> nu = {0.0, 0.0};
  std::vector<double> coef = {sqrt(2.0), sqrt(3.0), sqrt(4.0)};

  double integrand;
  /*std::vector<std::vector<double>> lambdas =  {{0.1543656365691809, 12.154365636569182},
											   {0.2483250102837163, 12.248325010283716},
											   {0.3422843839982518, 12.342284383998251},
											   {0.4362437577127873, 12.436243757712788},
											   {0.5302031314273228, 12.530203131427323},
											   {0.6241625051418582, 12.624162505141857},
											   {0.7181218788563937, 12.718121878856394},
											   {0.8120812525709291, 12.812081252570929},
											   {0.9060406262854646, 12.906040626285463},
											   {1.0, 13.0}};*/
  std::vector<std::vector<double>> lambdas =  {{mu[0], 11.0},
											   {mu[1], 11.0},
											   {mu[2], 11.0},
											   {mu[3], 11.0},
											   {mu[4], 11.0},
											   {mu[5], 11.0},
											   {mu[6], 11.0},
											   {mu[7], 11.0},
											   {mu[8], 11.0},
											   {mu[9], 12.0}};
  
  std::vector<std::vector<double>> quad_params;
  std::vector<double> nodes, weights;

  for(int k=0; k<mu.size(); k++){
    // Initialise dummy variable to hold the value of the approximated integral
    double dummy = 0.0;
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
  		integrand = gsl_sf_bessel_Jnu(mu[k], coef[0]*nodes[j])*gsl_sf_bessel_Jnu(nu[0], coef[1]*nodes[j])*gsl_sf_bessel_Jnu(nu[1], coef[2]*nodes[j]);
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
    /*std::cout << "k=" << k+1 << ",   In = " << std::setprecision(std::numeric_limits<double>::max_digits10) << In[k] 
    					               << ",   Tn = " << std::setprecision(std::numeric_limits<float>::max_digits10) << Tn[k] << " seconds\n";*/
	std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << In[k] << "\n";  
  }

  return 0;
}