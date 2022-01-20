//---------------------------------------------------------------------------------------
// File:      src/MonMap.cpp
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
//       FUNCTION: lambda_max = computeLambdaMax(lambda_min, user_n)
//                
//          INPUT: - lambda_min = minimum (and only) exponent in the input sequence
//                 - user_n = desired number of (quadrature) nodes provided by the user
//
//         OUTPUT: - lambda_max = additional exponent of the new binomial
//
//    DESCRIPTION: the monomial transformation quadrature rule processes the G-L nodes & 
//                 weights for those polynomials characterised by an arbitrarly large
//                 gap between the terms of minimum and maximum degree. There might be
//                 cases however in which the user wants to integrate singular 
//                 monomials; in those cases the library will require an additional, 
//                 non-constant, term to be added to the monomial. It does so by either 
//                 allowing the user to manually input the exponent of the additional
//                 term from the CLI or to specify the number of nodes to use in
//                 its application. In this last instance the following function 
//                 automatically generates the resulting exponent of the additional
//                 term (lambda_max).
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
//    DESCRIPTION: once the exponents in the terms with minimum and maximum degree in
//                 the user-input polynomial have been determined, this method 
//                 solves equation 1.18, which is a 7-th degree polynomial in
//                 n derived by a linear regression of n in beta_min/beta_max. 
//                 It then extracts the only real root whose integer floor will then be 
//                 the minimum possible number of (quadrature) nodes to be used in the
//                 new quadrature formula to achieve double precision (the solver itself
//                 is a class' method implemented in GSL-GNU Scientific Library).
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
//    DESCRIPTION: the order of the monomial transformation that characterises the new
//                 samples of the quadrature rule obtained starting from those of the
//                 G-L quadrature formula is computed as a linear interpolation
//                 beetween the minimum and maximum bound outlined in the inequality
//                 1.17 of the doc/UserManual.pdf.
//
/////////////////////////////////////////////////////////////////////////////////////////

double computeMapOrder(const std::vector<float128>& lambdas, const std::vector<float128>& betas)
{
  // Compute r as the linear interpolation between r_min and r_max
  float128 transf_order = ((static_cast<float128>(1.0) + betas[0])/(static_cast<float128>(1.0) + lambdas[0])+
                         (static_cast<float128>(1.0) + betas[1])/(static_cast<float128>(1.0) + lambdas[1]))/static_cast<float128>(2.0);
  
  // Print on-screen information
  std::cout << std::setprecision(std::numeric_limits<float>::max_digits10)
            << "\n ** Transformation order = "
            << transf_order << " **"
            << std::endl;

  // Output the results on text file
  std::string mkdir = "mkdir -p ";
  std::string output_dir = "output";
  std::string results_dir = mkdir + output_dir;
  system(results_dir.c_str());

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

  return static_cast<double>(transf_order);
}


/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: [{new_x, new_w, old_x, old_w}] = computeParamsGl(r, n_min)
//                
//          INPUT: - r = output of function 'computeMapOrder'
//                 - n_min = output of function 'streamMonMapData'
//
//         OUTPUT: - {new_x, new_w} = new set of G-L quadrature nodes (x) and
//                                    weights (w) following the monomial map
//                 - {old_x, old_w} = classical set of G-L quadrature nodes (x) and
//                                    weights (w) following the affine map
//
//    DESCRIPTION: once the transformation order is available, the monomial map itself
//                 is applied to the G-L nodes and weights that have been previously
//                 mapped from [-1,1] to [0,1] via an affine (linear) map. The 
//                 complete set of new and old nodes and weights (referred to as
//                 quadrature_parameters) are collected in a tuple and outputted by
//                 this method.
//
/////////////////////////////////////////////////////////////////////////////////////////

std::tuple<std::vector<float128>, std::vector<float128>, std::vector<float128>, std::vector<float128>> computeParamsGl(const double& r, const int& n_min)
{
  std::ofstream results;
  results.open("output/Results.txt", std::ios_base::app);

  results << "\n\n OUTPUTS:"
          << "\n        num. of nodes = "
          << n_min;

  results.close();

  std::vector<float128> new_nodes, new_weights, old_nodes, old_weights;

  // Compute affine map [a,b] -> [-1,1] parameters
  float128 a = 0.0, b = 1.0;
  float128 alpha = 0.5*(b - a), beta = 0.5*(a + b);
  float128 jacobian = alpha;

  // Derive new and old G-L nodes
  std::ifstream nodes_file;
  nodes_file.open("../data/TabulatedGlNodes.csv");

  std::string line_nodes, column_nodes;
  float128 affine_node;
  int num_nodes;

  // Loop through the rows of TabulatedGlNodes.csv
  while(std::getline(nodes_file, line_nodes))
  {
    std::stringstream column_nodes_string(line_nodes);
    std::vector<std::string> row_nodes;

    // Loop through the columns of TabulatedGlNodes.csv
    while(std::getline(column_nodes_string, column_nodes, ','))
    {
      row_nodes.push_back(column_nodes);
    }

    // Extract n as the first column of TabulatedGlNodes.csv
    num_nodes = stoi(row_nodes[0]);

    if(n_min == num_nodes)
    {
      for(int k = 1; k <= num_nodes; k++)
      {
        // Map the original G-L node from [-1,1] to [0,1] through affine transformation
        affine_node = alpha*(static_cast<float128>(row_nodes[k].substr(1, row_nodes[k].size() - 2))) + beta;
        old_nodes.push_back(affine_node);
        // Compute new G-L node using the monomial transformation
        new_nodes.push_back(pow(affine_node, static_cast<float128>(r)));
      }

      break;
    }
  }
  nodes_file.close();

  // Derive new and old G-L weights
  std::ifstream weights_file;
  weights_file.open("../data/TabulatedGlWeights.csv");

  std::string line_weights, column_weights;
  float128 affine_weight;
  int num_weights;

  // Loop through the rows of TabulatedGlWeights.csv
  while(std::getline(weights_file, line_weights))
  {
    std::stringstream column_weights_string(line_weights);
    std::vector<std::string> row_weights;

    // Loop through the columns of TabulatedGlWeights.csv
    while(std::getline(column_weights_string, column_weights, ','))
    {
      row_weights.push_back(column_weights);
    }

    // Extract n as the first column of TabulatedGlNodes.csv
    num_weights = stoi(row_weights[0].substr(1, row_weights[0].size() - 2));

    if(n_min == num_weights)
    {

      for(int k = 1; k <= num_nodes; k++)
      {
        // Map the original G-L weight from [-1,1] in [0,1] through affine transformation
        affine_weight = jacobian*(static_cast<float128>(row_weights[k].substr(1, row_weights[k].size() - 2)));
        old_weights.push_back(affine_weight);
        // Compute new G-L weight using the monomial transformation
        new_weights.push_back(static_cast<float128>(r)*pow(old_nodes.at(k-1), static_cast<float128>(r-1))*affine_weight);
      }

      break;
    }
  }
  weights_file.close();

  return std::make_tuple(new_nodes, new_weights, old_nodes, old_weights);
}


/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: In = computeQuadGl(x, w, muntz_sequence, coeff_sequence)
//                
//          INPUT: - x = output of function 'computeParamsGl' optimised by 'optimiseData'
//                 - w = output of function 'computeParamsGl' optimised by 'optimiseData'
//                 - muntz_sequence = sequence of real exponents of the polynomial
//                 - coeff_sequence = sequence of real coefficients of the polynomial
//
//         OUTPUT: - In = value of the numerical approximated integral for the user input
//
//    DESCRIPTION: every interpolatory quadrature rule approximates the definite integral
//                 by means of a weighted sum of the kernel's values on specific points
//                 along the integration interval (i.e. nodes); an interpolatory
//                 Gaussian quadrature formula is a quadrature rule whose nodes 
//                 corresponds to the roots of a polynomial that is orthogonal in the 
//                 integration interval to the weight function of the kernel. A G-L 
//                 quadrature formula with n+1 nodes is a Gaussian formula for which the
//                 nodes corresponds to the roots of the Legendre n-th degree polynomial
//                 that is orthogonal to the weight function w(x) = 1 in [-1, 1]. This 
//                 routine implements the computation of such weighted sum in a general 
//                 fashion, i.e. regardless of the quadrature rules, as it only requires
//                 the specified samples. It is thus re-used multiple times in other
//                 modules to compute both the classical G-L and the monomial transfor-
//                 mation quadrature rules.
//
/////////////////////////////////////////////////////////////////////////////////////////

template<typename type1, typename type2>
float128 computeQuadGl(const std::vector<type1>& nodes,const std::vector<type1>& weights, std::vector<type2>& muntz_sequence, std::vector<type2>& coeff_sequence)
{
  // The G-L quadrature retains the same precision as the input nodes and weights with which it is computed
  float128 In = static_cast<float128>(0.0);
  // Loop over the monomial terms of the polynomial
  for(int j=0; j < muntz_sequence.size(); j++)
  {
    // Vector containing the single terms of (x_j)^k for each node j=0,...,n_min
    std::vector<float128> f_values;
    // Loop over G-L nodes
    for(int k = 0; k < nodes.size(); k++)
    {
      // Compute the single term f(x_j) = (x_j)^lambda[k] and store it in the f_values vector
      f_values.push_back(pow(static_cast<float128>(nodes[k]), static_cast<float128>(muntz_sequence[j])));
    }
    In += static_cast<float128>(coeff_sequence[j])*doubleDotProduct(f_values, weights);
  }
  return In;
}
template float128 computeQuadGl<float128, float128>(const std::vector<float128>& nodes, const std::vector<float128>& weights, std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence); // Template mock instantiation for non-inline function
template float128 computeQuadGl<float128, double>(const std::vector<float128>& nodes, const std::vector<float128>& weights, std::vector<double>& muntz_sequence, std::vector<double>& coeff_sequence); // Template mock instantiation for non-inline function
template float128 computeQuadGl<double, float128>(const std::vector<double>& nodes, const std::vector<double>& weights, std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence); // Template mock instantiation for non-inline function
template float128 computeQuadGl<double, double>(const std::vector<double>& nodes, const std::vector<double>& weights, std::vector<double>& muntz_sequence, std::vector<double>& coeff_sequence); // Template mock instantiation for non-inline function


/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: En = computeExactError({post_map_quadrature, pre_map_quadrature}, 
//                                       muntz_sequence, coeff_sequence, print_flag)
//                
//          INPUT: - {post_map_quadrature, pre_map_quadrature} = output of function
//                                                               'computeQuadGl'
//                 - muntz_sequence = sequence of real exponents of the polynomial
//                 - coeff_sequence = sequence of real coefficients of the polynomial
//                 - print_flag = boolean parameter to tell the routine wheter or not
//                                to print the numerical value of the primitive
//
//         OUTPUT: - En = relative error of the specified rule computed with the new
//                        nodes and weights
//
//    DESCRIPTION: let I_n be the numerical integral calculated using the G-L quadrature
//                 rule (outputed by function 'computeQuadGl') and I_ex be the exact
//                 (analytic) value of such integral then the relative error of the 
//                 quadrature can be computed a-posteriori as R_n = |I_ex - I_n|/|I_ex|.
//                 Such computation is implemented in this routine with the exact 
//                 integral being as precise as the user-input polynomial is.
//
/////////////////////////////////////////////////////////////////////////////////////////

template<typename type>
float128 computeExactError(const float128& In, std::vector<type>& muntz_sequence, std::vector<type>& coeff_sequence, bool& print_primitive)
{
  // The exact integral's precision is defined by the precision of the original G-L parameters and the coefficients and exponents of the input polynomial which is at best float50
  float128 I = static_cast<float128>(0.0);
  // Computing the definite integral of each monomial in [0,1]
  for(int k=0; k < muntz_sequence.size(); k++)
  {
    // The constant term multiplies the difference (b^{k+1} - a^{k+1}) = 1 (since a=0 and b=1) which is thus omitted
    I += static_cast<float128>(coeff_sequence[k])/(static_cast<float128>(muntz_sequence[k]) + static_cast<float128>(1.0));
  }
  // Compute and return the exact a-posteriori remainder of the quadrature formula
  if(print_primitive)
  {
     std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
               << " ** I(p(x))   = "
               << I;
  }
  return fabs(I - In)/fabs(I);
}
template float128 computeExactError(const float128& In, std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence, bool& print_primitive);
template float128 computeExactError(const float128& In, std::vector<double>& muntz_sequence, std::vector<double>& coeff_sequence, bool& print_primitive);