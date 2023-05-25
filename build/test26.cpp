#include <gsl/gsl_sf_bessel.h>
#include <chrono>
#include "mtqr.h"

// clear && rm outputs.txt && g++ -w test26.cpp -o test26 -I../include -L. -lmtqr -lm -lquadmath -lgsl -lgslcblas && ./test26

int main(int argc, char** argv){
  
  double a = 1.0/2.0;
  double b = -1.0/2.0;

  std::vector<double> In;
  std::vector<double> Tn;
  
  std::vector<double> nu = {a, b};
  std::vector<double> coef = {0.1, 10.0};

  double integrand;
  std::vector<double> params = {0.001, 
                                0.006210526315789474, 
                                0.01142105263157895, 
                                0.01663157894736842, 
                                0.0218421052631579, 
                                0.02705263157894737, 
                                0.03226315789473684, 
                                0.03747368421052632, 
                                0.04268421052631579,
                                0.04789473684210526, 
                                0.05310526315789475, 
                                0.05831578947368422,
                                0.06352631578947368, 
                                0.06873684210526317, 
                                0.07394736842105264,
                                0.07915789473684211, 
                                0.08436842105263159, 
                                0.08957894736842106,
                                0.09478947368421053, 
                                0.1};
  std::vector<double> lambdas = {-3.0/4.0, 52.0};

  double m = b - a + 1.0;

  // Apply the monomial transformation and obtain the new, optimised quadrature rule's nodes and weights
  auto start = std::chrono::high_resolution_clock::now();
  std::vector<std::vector<double>> quad_params = mtqr(lambdas[0], lambdas[1]);
  auto stop = std::chrono::high_resolution_clock::now();

  std::vector<double> nodes = quad_params.at(0);
  std::vector<double> weights = quad_params.at(1);

  for(int k=0; k<params.size(); k++){
    // Initialise dummy variable to hold the value of the approximated integral
    double dummy = 0.0;
    double u = params[k];

  	// Debug
  	// std::cout << std::setprecision(std::numeric_limits<double>::max_digits10) << "k=" << k+1 << ", nu_1=" << nu_1 << ", nu_2=" << nu_2 << ", a_1=" << coeff_1 << ", a_2=" << coeff_2 << "\n";
  	
    // Compute the definite integral numerically using the MTQR with new nodes and weights
  	for(int j=0; j<nodes.size(); j++){
  		// Compute the value of the Bessel function of (fractional order) nu=nu_min when evaluated @ the j-th node of the loop
  		integrand = (pow(nodes[j], m)/(pow(nodes[j], 2.0) + pow(u, 2.0)))*gsl_sf_bessel_Jnu(a, coef[0]*nodes[j])*gsl_sf_bessel_Jnu(b, coef[1]*nodes[j]);
  		// Update the quadrature approximation with the current contribution
  		dummy += weights[j]*integrand;
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