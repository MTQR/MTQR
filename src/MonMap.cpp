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
//       FUNCTION: E_n = computeEstimate(lambda, n, envelope)
//                
//          INPUT: - lambda = see (18) in [1]
//                 - n = see (18) in [1]
//                 - envelope = string flag that specifies whether the estimate must be 
//                              enveloped (to avoid finite arithmetic infinities) or not
//
//         OUTPUT: - E_n = exact asymptotic estimate of the G-L quadrature error 
//                   computed using (13) and (18) in [1]
//
//    DESCRIPTION: in order for the library to generate tabulated values for 
//                 beta_min/beta_max as functions of n the exact estimates of the G-L 
//                 quadrature error must be computed with high-precision. Such an 
//                 a-priori estimate is generally known as results of complex analysis
//                 (see [2]) however in [1] a more accurate form has been devised in
//                 formulae (13) and (18) which are implemented in this routine to 
//                 compute the aformentioned error estimate.
//
//      REFERENCE: [1] = Lombardi Guido - Design of quadrature rules for Müntz and 
//                                        Müntz-logarithmic polynomials using monomial
//                                        transformation,
//                                        Int. J. Numer. Meth. Engng., 80: 1687-1717,
//                                        https://doi.org/10.1002/nme.2684.
//                 [2] = Donaldson J.D., Elliott D. - A unified approach to quadrature
//                                        rules with asymptotic estimates of their
//                                        remainders,
//                                        SIAM Journal on Numerical Analysis 1972;
//                                        9(4):573–602,
//                                        https://doi.org/10.1137/0709051
//
/////////////////////////////////////////////////////////////////////////////////////////

float1k computeErrorEstimate(const float128& input_lambda, const int& num_nodes, const std::string& envelope)
{

  float1k lambda = static_cast<float1k>(input_lambda);

  float1k common_factor = 2*(boost::math::powm1(2,-(1+lambda))+1);
  float1k denominator = 2*num_nodes + lambda;

  float1k b1 = (boost::math::tgamma(2*lambda)*boost::math::tgamma(2*num_nodes - lambda))/boost::math::tgamma(2*lambda + 2*num_nodes - lambda);
  float1k b2 = (boost::math::tgamma(2*lambda)*boost::math::tgamma(2 + 2*num_nodes - lambda))/boost::math::tgamma(2*lambda + 2 + 2*num_nodes - lambda);

  float1k exact = common_factor*(boost::math::powm1(2,-lambda)+1)*lambda*((b1/denominator) - b2/(2 + denominator));
  
  if(envelope.compare("yes") == 0)
  {
    if(lambda + 1/2 - 2*num_nodes >= 0)
    {
      return fabs(exact*boost::math::sin_pi(lambda));
    }
    else
    {
      return fabs(exact);
    }
  }
  else
  {
    return fabs(exact*boost::math::sin_pi(lambda));
  }
}


/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: lambda_max = computeLambda(lambda_min, user_n)
//                
//          INPUT: - lambda_min = minimum exponent in the "muntz_sequence" input of 
//                   function 'getInputData'
//                 - user_n = desired number of (quadrature) nodes defined by the user
//
//         OUTPUT: - lambda_max = computed additional exponent of the polynomial
//
//    DESCRIPTION: the monomial quadrature rule is a pre-processing of the G-L nodes & 
//                 weights for those polynomials characterised by an arbitrarly large
//                 gap between the terms of minimum and maximum degree. There might be
//                 cases however in which the user wants to integrate singular 
//                 monomials; in those cases the library will require an additional, 
//                 non-constant, term to be added to the monomial. It does so by either 
//                 allowing the user to either manually input the exponent of the 
//                 additional term from the CLI (see lines 180~184 in the src/DatIo.cpp
//                 file) or to specify the maximum number of quadrature nodes to use in
//                 its application. In this last instance the the following function 
//                 automatically generates the resulting exponent of the additional
//                 term (lambda_max).
//
/////////////////////////////////////////////////////////////////////////////////////////

float128 computeLambdaMax(float128& lambda_min, int num_nodes)
{
  auto data = streamMonMapData(num_nodes);
  
  int n_max = std::get<0>(data);
  std::vector<double> betas = std::get<1>(data);

  double r_min = (1 + betas[0])/(1 + static_cast<double>(lambda_min));
  float128 lambda_max = (1 + betas[1] - r_min)/r_min;
  
  return lambda_max;
}


/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: n = computeNumNodes(lambda_min, lambda_max)
//                
//          INPUT: - lambda_min = minimum exponent in the input "muntz_sequence" 
//                   (strictly greater than -1)
//                 - lambda_max = maximum exponent in the input "muntz_sequence"
//
//         OUTPUT: - n = number of (quadrature) nodes computed as solution of equation
//                      (62) in [1]
//
//    DESCRIPTION: once the exponents in the terms with minimum and maximum degree in
//                 the user-input polynomial have been determined, this method 
//                 implements formula (62) in [1], which is a 7-th degree polynomial in
//                 n, and extract its only real root whose integer floor will then be 
//                 the minimum possible number of (quadrature) nodes to be used in G-L
//                 formula to achieve double precision (the polynomial solver itself
//                 is a class method implemented in GSL-GNU Scientific Library).
//
//      REFERENCE: [1] = Lombardi Guido - Design of quadrature rules for Müntz and 
//                                        Müntz-logarithmic polynomials using monomial
//                                        transformation,
//                                        Int. J. Numer. Meth. Engng., 80: 1687-1717,
//                                        https://doi.org/10.1002/nme.2684.
//
/////////////////////////////////////////////////////////////////////////////////////////

int computeNumNodes(const float128& lambda_min, const float128& lambda_max)
{
  double l_max = static_cast<double>(lambda_max), l_min = static_cast<double>(lambda_min); 
  constexpr double c0 = -0.0040693, c1 = 0.00041296, d0 = 7.8147, d2 = 0.10123;

  double c2 = d2 + d2*l_min;
  double c3 = d0 + d0*l_min + l_min - l_max;
  double c4 = pow(1+l_max,3.0) + 1;

  double coeff7 = static_cast<double>(c1*(pow(c2,3)));
  double coeff6 = static_cast<double>(c0*(pow(c2,3)));
  double coeff5 = static_cast<double>(3*c1*(pow(c2,2))*c3);
  double coeff4 = static_cast<double>(3*c0*(pow(c2,2))*c3);
  double coeff3 = static_cast<double>(3*c1*c2*(pow(c3,2)));
  double coeff2 = static_cast<double>(3*c0*c2*(pow(c3,2)));
  double coeff1 = static_cast<double>(c1*(pow(c3,3)));
  double coeff0 = static_cast<double>(c1*(pow(c3,3))-c4);

  double coeff[8] = {coeff0, coeff1, coeff2, coeff3, coeff4, coeff5, coeff6, coeff7};
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
//       FUNCTION: r = computeOrder({lambda_min, lambda_max}, {beta_min, beta_max})
//                
//          INPUT: - {lambda_min, lambda_max} = output of function 'getInputData'
//                 - {beta_min, beta_max} = output of function 'retrieveMonData'
//
//         OUTPUT: - r = real value of the transformation order of the monomial map
//
//    DESCRIPTION: the order of the monomial transformation for the new nodes and weight
//                 of the G-L quadrature formula is computed as a linear interpolation
//                 beetween r_min and r_max reported in (63) of [1].
//
//      REFERENCE: [1] = Lombardi Guido - Design of quadrature rules for Müntz and 
//                                        Müntz-logarithmic polynomials using monomial
//                                        transformation,
//                                        Int. J. Numer. Meth. Engng., 80: 1687-1717,
//                                        https://doi.org/10.1002/nme.2684.
//
/////////////////////////////////////////////////////////////////////////////////////////


double computeMapOrder(const std::vector<float128>& lambdas, const std::vector<double>& betas)
{
  double transf_order = ((1 + betas[0])/(1 + static_cast<double>(lambdas[0]))+
                         (1 + betas[1])/(1 + static_cast<double>(lambdas[1])))/2;
  
  std::cout << std::setprecision(std::numeric_limits<float>::max_digits10)
            << "\n ** Transformation order = "
            << transf_order << " **"
            << std::endl;

  std::string mkdir = "mkdir -p ";
  std::string output_dir = "output";
  std::string results_dir = mkdir + output_dir;
  system(results_dir.c_str());

  std::string results_file = output_dir + "/Results.txt";
  std::ofstream results;
  results.open(results_file);

  results << std::setprecision(std::numeric_limits<float128>::max_digits10)
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

  return transf_order;
}


/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: [J, {new_x, new_w, old_x, old_w}] = computeParams(r, n_min)
//                
//          INPUT: - r = output of function 'computeOrder'
//                 - n_min = output of function 'retrieveMonData'
//
//         OUTPUT: - J = jacobian of the affine map phi: [a,b] -> [-1,1] of the G-L 
//                       quadrature formula
//                 - {new_x, new_w} = new set of G-L quadrature nodes (x) and
//                                    weights (w) following the monomial map
//                 - {old_x, old_w} = classical set of G-L quadrature nodes (x) and
//                                    weights (w) in [-1, 1]
//
//    DESCRIPTION: once the transformation order is available, the monomial map itself
//                 is constructed according to (55) in [1]; of course, prior to the 
//                 monomial map is applied, the G-L nodes and weights have to be mapped
//                 from the user-input interval I = [a,b] to [-1, 1]. This is 
//                 implemented in this routine which also carries the jacobian of the
//                 affine map among the outputs since it is the multiplication
//                 coefficient of the integral itself (thus needed by the function
//                 'computeQuadGl'). Furthermorethe classic G-L nodes and weights are
//                 outputted as well so that comparisons ca be made between the
//                 traditional G-L quadrature integral and the monomial rule quadrature.
//
/////////////////////////////////////////////////////////////////////////////////////////

std::tuple<std::vector<float50>, std::vector<float50>, std::vector<float50>, std::vector<float50>> computeParamsGl(const double& r, const int& n_min)
{
  std::ofstream results;
  results.open("output/Results.txt", std::ios_base::app);

  results << "\n\n OUTPUTS:"
          << "\n        num. of nodes = "
          << n_min;

  results.close();

  std::vector<float50> new_nodes, new_weights, old_nodes, old_weights;

  // Compute affine map [a,b] -> [-1,1] parameters
  constexpr double a = 0, b = 1;
  double alpha = 0.5*(b - a), beta = 0.5*(a + b);
  double jacobian = alpha;

  // Derive new and old G-L nodes
  std::ifstream nodes_file;
  nodes_file.open("../data/TabulatedGlNodes.csv");

  std::string line_nodes, column_nodes;
  float50 affine_node;
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
        affine_node = alpha*(static_cast<float50>(static_cast<float128>(static_cast<float50>(row_nodes[k].substr(1, row_nodes[k].size() - 2))))) + beta;
        old_nodes.push_back(affine_node);
        // Compute new G-L node using the monomial transformation
        new_nodes.push_back(static_cast<float50>(pow(affine_node,r)));
      }

      break;
    }
  }
  nodes_file.close();

  // Derive new and old G-L weights
  std::ifstream weights_file;
  weights_file.open("../data/TabulatedGlWeights.csv");

  std::string line_weights, column_weights;
  float50 affine_weight;
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
        affine_weight = jacobian*(static_cast<float50>(static_cast<float128>(static_cast<float50>(row_weights[k].substr(1, row_weights[k].size() - 2)))));
        old_weights.push_back(affine_weight);
        // Compute new G-L weight using the monomial transformation
        new_weights.push_back(static_cast<float50>(r*(pow(old_nodes.at(k-1), r-1))*affine_weight));
      }

      break;
    }
  }
  weights_file.close();

  return std::make_tuple(new_nodes, new_weights, old_nodes, old_weights);
}


/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: quadrature = computeQuadGl(x, w, muntz_sequence, coeff_sequence)
//                
//          INPUT: - x = G-L quadrature nodes
//                 - w = G-L quadrature weights
//                 - muntz_sequence = sequence of real exponents of the polynomial
//                 - coeff_sequence = sequence of real coefficients of the polynomial
//
//         OUTPUT: - quadrature = G-L quadrature of user-input polynomial
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
//                 routine implements the G-L quadrature formula provided classical and
//                 new G-L nodes and weights.
//
/////////////////////////////////////////////////////////////////////////////////////////

template<typename type>
type computeQuadGl(const std::vector<type>& nodes,const std::vector<type>& weights, std::vector<type>& muntz_sequence, std::vector<type>& coeff_sequence)
{
  // The G-L quadrature retains the same precision as the input nodes and weights with which it is computed
  type In = 0;
  // Loop over the monomial terms of the polynomial
  for(int j=0; j < muntz_sequence.size(); j++)
  {
    // Vector containing the single terms of (x_j)^k for each node j=0,...,n_min
    std::vector<type> f_values;
    // Loop over G-L nodes
    for(int k = 0; k < nodes.size(); k++)
    {
      // Compute the single term f(x_j) = (x_j)^lambda[k] and store it in the vector
      f_values.push_back(pow(nodes[k], muntz_sequence[j]));
    }
    In += coeff_sequence[j]*orderedInnerProduct(f_values, weights);
  }
  return In;
}
template float50 computeQuadGl<float50>(const std::vector<float50>& nodes, const std::vector<float50>& weights, std::vector<float50>& muntz_sequence, std::vector<float50>& coeff_sequence); // Template mock instantiation for non-inline function
template float128 computeQuadGl<float128>(const std::vector<float128>& nodes, const std::vector<float128>& weights, std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence); // Template mock instantiation for non-inline function
template double computeQuadGl<double>(const std::vector<double>& nodes, const std::vector<double>& weights, std::vector<double>& muntz_sequence, std::vector<double>& coeff_sequence); // Template mock instantiation for non-inline function


/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: error = computeError({post_map_quadrature, pre_map_quadrature}, 
//                                       muntz_sequence, coeff_sequence, I)
//                
//          INPUT: - {post_map_quadrature, pre_map_quadrature} = output of function
//                                                               'computeQuadGl'
//                 - muntz_sequence = sequence of real exponents of the polynomial
//                 - coeff_sequence = sequence of real coefficients of the polynomial
//                 - I = [a,b] = interval of integration of the user-input polynomial
//
//         OUTPUT: - error = relative error of the G-L quadrature computed with the new
//                           nodes and weights
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
type computeExactError(const type& In, std::vector<type>& muntz_sequence, std::vector<type>& coeff_sequence)
{
  // The exact integral's precision is defined by the precision of the original G-L parameters and the coefficients and exponents of the input polynomial which is at best float50
  type I = 0;
  // Computing the definite integral of each monomial in [0,1]
  for(int k=0; k < muntz_sequence.size(); k++)
  {
    // The constant term multiplies the difference (b^{k+1} - a^{k+1}) = 1 (since a=0 and b=1) which is thus omitted
    I += coeff_sequence[k]*(1/(muntz_sequence[k]+1));
  }
  // Compute and return the exact a-posteriori remainder of the quadrature formula
  std::cout << "\nI(f) = " << std::setprecision(std::numeric_limits<float50>::max_digits10) << I << std::endl;
  return fabs(I - In)/fabs(I);
}
template float50 computeExactError<float50>(const float50& In, std::vector<float50>& muntz_sequence, std::vector<float50>& coeff_sequence); // Template mock instantiation for non-inline function
template float128 computeExactError<float128>(const float128& In, std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence); // Template mock instantiation for non-inline function
template double computeExactError <double>(const double& In, std::vector<double>& muntz_sequence, std::vector<double>& coeff_sequence); // Template mock instantiation for non-inline function