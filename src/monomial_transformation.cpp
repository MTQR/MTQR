//---------------------------------------------------------------------------------------
// File:      src/monomial_transformation.cpp
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

#include "mtqr.h"

#include "GL_NODES.h"   // Include Gauss-Legendre nodes in (-1,1) with 50 decimal digits
#include "GL_WEIGHTS.h" // Include Gauss-Legendre weights in (0,1) with 50 decimal digits

/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: lambda_max = computeLambdaMax(lambda_min, user_n)
//                
//          INPUT: - lambda_min = minimum (and only) exponent in the input sequence
//                 - user_n = desired number of (quadrature) nodes provided by the user
//
//         OUTPUT: - lambda_max = additional exponent of the new binomial
//
//    DESCRIPTION: there might be cases in which the user wants to integrate generalised 
//                 monomials rather than the alike polynomials. In those cases the 
//                 monomial transformation quadrature rule requires additional input to 
//                 work properly. The library handles such case by either allowing the
//                 user to manually input the exponent of the additional term (from the 
//                 CLI) or by letting the user specify the number of nodes to be used.
//                 In this last instance this method automatically computes the resulting
//                 maximum exponent of an additional term (lambda_max) that the monomial
//                 transformation quadrature rule can precisely integrate.
//
/////////////////////////////////////////////////////////////////////////////////////////

float128 computeLambdaMax(float128& lambda_min, int num_nodes)
{
  auto data = streamMonMapData(num_nodes);
  
  int n_max = std::get<0>(data);
  std::vector<float128> betas = std::get<1>(data);

  float128 r_min = (static_cast<float128>(1.0) + betas[0])/(static_cast<float128>(1.0) + static_cast<float128>(lambda_min));
  float128 lambda_max = (static_cast<float128>(1.0) + betas[1] - r_min)/r_min;
  
  return lambda_max;
}


/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: n = computeNumNodes(lambda_min, lambda_max)
//                
//          INPUT: - lambda_min = minimum exponent in the input "muntz_sequence" 
//                                (strictly greater than -1)
//                 - lambda_max = maximum exponent in the input "muntz_sequence"
//
//         OUTPUT: - n = number of (quadrature) nodes computed as the only real solution
//                       of equation 1.18 reported in the doc/UserManual.pdf
//
//    DESCRIPTION: this method solves equation (1.18) in doc/UserManual.pdf, which is a 
//                 7-th degree polynomial in n_min (the number of quadrature nodes).
//                 Of the 7 roots, it then extracts the only real one and then takes its
//                 ceil to return an integer value for n_min. The solver itself is a 
//                 class' method implemented in the GSL-GNU Scientific Library.
//
/////////////////////////////////////////////////////////////////////////////////////////

int computeNumNodes(const float128& lambda_min, const float128& lambda_max)
{
  constexpr double c0 = -0.0040693, c1 = 0.00041296, d0 = 7.8147, d2 = 0.10123;

  double c2 = d2 + d2*static_cast<double>(lambda_min);
  double c3 = d0 + d0*static_cast<double>(lambda_min) + static_cast<double>(lambda_min) - static_cast<double>(lambda_max);
  double c4 = pow(1.0 + static_cast<double>(lambda_max), 3.0);

  double coeff[8] = {c1*pow(c3,3.0) - c4,
                     c1*pow(c3,3.0),
                     3.0*c0*c2*pow(c3,2.0),
                     3.0*c1*c2*pow(c3,2.0),
                     3.0*c0*pow(c2,2.0)*c3,
                     3.0*c1*pow(c2,2.0)*c3,
                     c0*pow(c2,3.0),
                     c1*pow(c2,3.0)};
  double n[14];

  gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc (8);
  gsl_poly_complex_solve(coeff, 8, w, n);
  gsl_poly_complex_workspace_free(w);

  int n_min, flag = 0;

  for (int k = 0; k < 7; k++)
  {
    if(fabs(n[2*k+1]) < EPS)
    {
      if(n[2*k] > 0 && flag==0)
      {
        n_min = ceil(n[2*k]);
        flag = 1;
      }
    }
  }

  return n_min;
}


/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: r = computeMapOrder({lambda_min, lambda_max}, {beta_min, beta_max})
//                
//          INPUT: - {lambda_min, lambda_max} = output of function 'manageData'
//                 - {beta_min, beta_max} = output of function 'streamMonMapData'
//
//         OUTPUT: - r = value of the transformation order of the monomial map (gamma)
//
//    DESCRIPTION: this method computes the order of the monomial transformation (1.16)
//                 in doc/UserManual.pdf. It is derived as a linear interpolation of the
//                 minimum and maximum bound of inequality (1.17) of doc/UserManual.pdf.
//
/////////////////////////////////////////////////////////////////////////////////////////

double computeMapOrder(const std::vector<float128>& lambdas, const std::vector<float128>& betas)
{
  // Compute r as the linear interpolation between r_min and r_max
  float128 transf_order = ((static_cast<float128>(1.0) + betas[0])/(static_cast<float128>(1.0) + lambdas[0])+
                         (static_cast<float128>(1.0) + betas[1])/(static_cast<float128>(1.0) + lambdas[1]))/static_cast<float128>(2.0);
  extern bool loud_mode;
  if(loud_mode)
  {
    // Create output sub-directory
    std::string mkdir = "mkdir -p ";
    std::string output_dir = "output";
    std::string results_dir = mkdir + output_dir;
    system(results_dir.c_str());

    // Print on-screen information
    std::cout << std::setprecision(std::numeric_limits<float>::max_digits10)
              << "\n ** Transformation order = "
              << transf_order << " **"
              << std::endl;

    // Output the results on text file
    std::string results_file = output_dir + "/Results.txt";
    std::ofstream results;
    
    results.open(results_file);
    results << std::setprecision(std::numeric_limits<float>::max_digits10)
            << "\nINPUTS:"
            << "\n        lambda_min = "
            << lambdas[0]
            << "\n        lambda_max = "
            << lambdas[1];
    results << std::setprecision(std::numeric_limits<float>::max_digits10)
            << "\n\nMONOMIAL MAP:"
            << "\n        beta_min = "
            << betas[0]
            << "\n        beta_max = "
            << betas[1]
            << "\n        transformation order = "
            << transf_order;
    results.close();
  }
  else
  {
    std::ofstream Outputs;
    Outputs.open("outputs.txt", std::ios_base::app);
    Outputs << std::setprecision(std::numeric_limits<float128>::max_digits10)
            << "Monomial Transformation Quadrature Rule:"
            << "\n    lambda_min = "
            << lambdas[0]
            << "\n    lambda_max = "
            << lambdas[1];
  }

  return static_cast<double>(transf_order);
}


/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: [{new_x, new_w, old_x, old_w}] = computeQuadParams(r, n_min)
//                
//          INPUT: - r = output of function 'computeMapOrder'
//                 - n_min = output of function 'streamMonMapData'
//
//         OUTPUT: - {new_x, new_w} = new set of monomial transformation quadrature rule
//                                    nodes (x) and weights (w)
//                 - {old_x, old_w} = old set of classical G-L quadrature nodes (x) and
//                                    weights (w) following the affine map
//
//    DESCRIPTION: with the transformation order being set, the monomial map is applied
//                 to transform the G-L nodes and weights, previously mapped from (-1,1)
//                 to (0,1) via the affine (linear) map in (1.12) in doc/UserManual.pdf.  
//                 The complete set of new and old nodes and weights are collected in a 
//                 tuple and returned.
//
/////////////////////////////////////////////////////////////////////////////////////////

std::tuple<std::vector<float128>, std::vector<float128>, std::vector<float128>, std::vector<float128>> computeQuadParams(const double& r, const int& n_min)
{
  extern bool loud_mode;
  if(loud_mode)
  {
    std::ofstream results;
    results.open("output/Results.txt", std::ios_base::app);

    results << "\n\n OUTPUTS:"
            << "\n        num. of nodes = "
            << n_min;

    results.close();
  }

  std::vector<float128> new_nodes, new_weights, old_nodes, old_weights;

  // Compute affine map [a,b] -> [-1,1] parameters
  float128 a = 0.0, b = 1.0;
  float128 alpha = 0.5*(b - a), beta = 0.5*(a + b);
  float128 jacobian = alpha;

  std::string line_nodes, column_nodes;
  float128 affine_node;
  int num_nodes;

  int it = 0;
  while(true)
  {
    std::vector<std::string> _nodes = GL_NODES[it];
    num_nodes = stoi(_nodes[0]);

    if(n_min == num_nodes)
    { 
      for(int k = 1; k <= num_nodes; k++)
      {
        // Map the original G-L node from [-1,1] to [0,1] through affine transformation
        affine_node = alpha*(static_cast<float128>(_nodes[k])) + beta;
        old_nodes.push_back(affine_node);
        // Compute new node using the monomial transformation
        new_nodes.push_back(pow(affine_node, static_cast<float128>(r)));
      }

      break;
    }
    else {
        it++;
    }
  }

  std::string line_weights, column_weights;
  float128 affine_weight;
  int num_weights;

  it = 0;
  while(true)
  {
    std::vector<std::string> _weights = GL_WEIGHTS[it];
    num_weights = stoi(_weights[0]);

    if(n_min == num_weights)
    {

      for(int k = 1; k <= num_nodes; k++)
      {
        // Map the original G-L weight from [-1,1] in [0,1] through affine transformation
        affine_weight = jacobian*(static_cast<float128>(_weights[k]));
        old_weights.push_back(affine_weight);
        // Compute new weight using the monomial transformation
        new_weights.push_back(static_cast<float128>(r)*pow(old_nodes.at(k-1), static_cast<float128>(r-1))*affine_weight);
      }

      break;
    }
    else {
        it++;
    }
  }

  return std::make_tuple(new_nodes, new_weights, old_nodes, old_weights);
}

/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: In = computeQuadRule(x, w, muntz_sequence, coeff_sequence)
//                
//          INPUT: - x = output of 'computeQuadParams' as optimised by 'optimiseData'
//                 - w = output of 'computeQuadParams' as optimised by 'optimiseData'
//                 - muntz_sequence = sequence of real exponents of the polynomial
//                 - coeff_sequence = sequence of real coefficients of the polynomial
//
//         OUTPUT: - In = value of the numerical approximated integral for the user input
//
//    DESCRIPTION: the computation of the numerical integral using a quadrature formula 
//                 involves the sum of the weighted values that the integrand assumes on
//                 specified samples of the interval. This method uses such information
//                 to implement the weighted sum and returns it. Note that this method
//                 does not implement a specific quadrature rule, but rather a generic
//                 one. In the context of the library it is therefore used for both the
//                 monomial transformation and the classical G-L quadrature rules.
//
/////////////////////////////////////////////////////////////////////////////////////////

template<typename T, typename U>
float128 computeQuadRule(const std::vector<T>& nodes,const std::vector<T>& weights, std::vector<U>& muntz_sequence, std::vector<U>& coeff_sequence)
{
  // The quadrature rule retains the same precision as the input nodes and weights with which it is computed
  float128 In = static_cast<float128>(0.0);
  // Loop over the monomial terms of the polynomial
  for(int j=0; j < muntz_sequence.size(); j++)
  {
    // Vector containing the single terms of (x_j)^k for each node j=0,...,n_min
    std::vector<float128> f_values;
    // Loop over the input nodes
    for(int k = 0; k < nodes.size(); k++)
    {
      // Compute the single term f(x_j) = (x_j)^lambda[k] and store it in the f_values vector
      f_values.push_back(pow(static_cast<float128>(nodes[k]), static_cast<float128>(muntz_sequence[j])));
    }
    In += static_cast<float128>(coeff_sequence[j])*doubleDotProduct(f_values, weights);
  }
  return In;
}
template float128 computeQuadRule<float128, float128>(const std::vector<float128>& nodes, const std::vector<float128>& weights, std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence); // Template mock instantiation for non-inline function
template float128 computeQuadRule<float128, double>(const std::vector<float128>& nodes, const std::vector<float128>& weights, std::vector<double>& muntz_sequence, std::vector<double>& coeff_sequence); // Template mock instantiation for non-inline function
template float128 computeQuadRule<double, float128>(const std::vector<double>& nodes, const std::vector<double>& weights, std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence); // Template mock instantiation for non-inline function
template float128 computeQuadRule<double, double>(const std::vector<double>& nodes, const std::vector<double>& weights, std::vector<double>& muntz_sequence, std::vector<double>& coeff_sequence); // Template mock instantiation for non-inline function


/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: En = computeExactError({post_map_quadrature, pre_map_quadrature}, 
//                                       muntz_sequence, coeff_sequence, print_flag)
//                
//          INPUT: - {post_map_quadrature, pre_map_quadrature} = output of function
//                                                               'computeQuadRule'
//                 - muntz_sequence = sequence of real exponents of the polynomial
//                 - coeff_sequence = sequence of real coefficients of the polynomial
//                 - print_flag = boolean parameter to tell the routine wheter or not
//                                to print the numerical value of the primitive
//
//         OUTPUT: - En = relative error of the specified rule computed with the new
//                        nodes and weights
//
//    DESCRIPTION: let I_n be the numerical integral calculated using the above method
//                 computeQuadRule for a specific quadrature rule and let I be the exact
//                 (analytic) counterpart of such integral. The relative error of the 
//                 quadrature can be computed a-posteriori as R_n = |I - I_n|/|I| as in.
//                 (1.7) of doc/UserManual.pdf. This method is the implementation of such
//                 error calculation being as precise as the user-input polynomial is.
//
/////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
float128 computeExactError(const float128& In, std::vector<T>& muntz_sequence, std::vector<T>& coeff_sequence, bool& print_primitive)
{
  extern bool loud_mode;
  float128 I = static_cast<float128>(0.0);
  // Computing the definite integral of each monomial in [0,1]
  for(int k=0; k < muntz_sequence.size(); k++)
  {
    // The constant term multiplies the difference (b^{k+1} - a^{k+1}) = 1 (since a=0 and b=1) which is thus omitted
    I += static_cast<float128>(coeff_sequence[k])/(static_cast<float128>(muntz_sequence[k]) + static_cast<float128>(1.0));
  }
  // Compute and return the exact a-posteriori remainder of the quadrature formula
  if(print_primitive && loud_mode)
  {
     std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
               << " ** I(p(x))   = "
               << I;
  }
  return fabs(I - In)/fabs(I);
}
template float128 computeExactError(const float128& In, std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence, bool& print_primitive);
template float128 computeExactError(const float128& In, std::vector<double>& muntz_sequence, std::vector<double>& coeff_sequence, bool& print_primitive);