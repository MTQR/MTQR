//---------------------------------------------------------------------------------------
// File:      src/data_management.cpp
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

#include "BETAS.h" // Include source data relative to bounded values of minimum and maximum betas in high-precision

/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: [n, {lambda_min, lambda_max}] = manageData(muntz_sequence, 
//                                                              coeff_sequence)
//                
//        INPUT: - muntz_sequence = sequence of real exponents of the polynomial
//               - coeff_sequence = sequence of real coefficients of the polynomial
//
//       OUTPUT: - n = output of function 'computeNumNodes' or input by the user
//               - lambda_min = minimum exponent in the input "muntz_sequence"
//               - lambda_max = maximum exponent in the input "muntz_sequence"
//
//    DESCRIPTION: the user-input polynomial is provided to the library, via Main.cpp, by.
//                 two unique lists of coefficients and exponents. This method first 
//                 checks the proper form of those inputs. In particular:
//                    - the number of exponents and the number of coefficients coincide;
//                    - the input polynomial is at least a binomial (otherwise the 
//                      further user-input as reported in doc/UserManual.pdf);
//                    - lambda_min > -1 (otherwise the input sequence of exponents is 
//                      not a Muntz sequence and the polynomial is not L1);
//                 Following the checks the exponents' sequence is then sorted locally
//                 to identify and return lambda_min/lambda_max at global level alongside
//                 the associated minimum number of nodes computed n_min.
//
/////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
std::tuple<int, std::vector<float128>> manageData(std::vector<T>& muntz_sequence, std::vector<T>& coeff_sequence)
{
  int num_nodes, n_min;
  float128 lambda_max, lambda_min;
  bool compute_n_min = true;
  extern bool loud_mode;

  if(loud_mode)
  {
    // Print initial message and input polynomial
    std::cout << std::endl;
    std::cout << "    |―――――――――――――――――――――――――――――――――――――――――――――――――|\n"
              << "    |                    ** MTQR **                   |\n"
              << "    |  ** MONOMIAL TRANSFORMATION QUADRATURE RULE **  |\n"
              << "    |―――――――――――――――――――――――――――――――――――――――――――――――――|\n";

    if(muntz_sequence.size()==coeff_sequence.size())
    {
        std::cout << "\n\n Input polynomial p(x) = ";
      for(int k=0; k < muntz_sequence.size(); k++)
      {
        if(coeff_sequence[k]>0)
        {
          std::cout << " +";
        }
        else
        {
          std::cout << " ";
        }

        if(coeff_sequence[k]==1)
        {
          std::cout << std::setprecision(std::numeric_limits<float>::max_digits10)
                    << "x^("
                    << muntz_sequence[k] << ") ";
        }
        else
        {
          std::cout << std::setprecision(std::numeric_limits<float>::max_digits10)
                    << coeff_sequence[k] << "*x^(" << muntz_sequence[k] << ")";
        }
      }
      std::cout << std::endl;
    }
    else
    {
      std::cout << "\n   ** ERROR ** The number of exponents doesn't match the number of coefficients in the input polynomial\n";
      exit(1);
    }

    // Check compliance of the input Muntz sequence
    if(muntz_sequence.size()==1)
    {
      if(muntz_sequence[0] <= -1.0)
      {
        std::cout << "\n   ** ERROR ** Lambda_min has to be (strictly) greater than -1 to be in a Muntz sequence\n";
        exit(1);
      }

      float128 additional_lambda;

      std::string input;
      std::cout << "\n ** WARNING ** Your input is a monomial of non-integer degree."
                << "\n               MTQR needs a binomial for double-precision quadrature."
                << "\n               How do you proceed? ['nodes' for n_min ~ 'lambda' for lambda_max]"
                << "\n               Input: ";
      std::cin >> input;

      if(input.compare("nodes") == 0)
      {
        std::cout << "\n\nPlease specify the desired number of quadrature nodes (number must be even): ";
        std::cin >> num_nodes;

        if(num_nodes % 2 == 0)
        {
          compute_n_min = false;

          float128 local_lambda_min = static_cast<float128>(muntz_sequence[0]);
          additional_lambda = computeLambdaMax(local_lambda_min, num_nodes);
        }
        else
        {
          std::cout << "\n   ** ERROR ** The number of nodes (n_min) must be even\n";
          exit(1);
        }
      }
      else if(input.compare("lambda") == 0)
      {
        std::cout << "\n\nPlease enter the exponent value with '.' separating the decimal digits from the integer part [the more decimal digits the better the precision is]: ";
        std::cin >> additional_lambda;
      }
      else
      {
        std::cout << "\n   ** ERROR ** Your Input must be either 'nodes' or 'lambda'"
                  << "\n               Please refer to Section 2.4 of doc/UserManual.pdf\n";
        exit(1);
      }

      muntz_sequence.push_back(static_cast<T>(additional_lambda));
      coeff_sequence.push_back(1.0);
    }
  }

  // Sort the input Muntz sequence to extract lambda_min and lambda_max
  std::vector<T> loc_muntz_seq = muntz_sequence;
  sort(loc_muntz_seq.begin(), loc_muntz_seq.end());
  if(loc_muntz_seq[0] <= -1.0)
  {
    std::cout << "\n   ** ERROR ** Lambda_min has to be (strictly) greater than -1 to be in a Muntz sequence\n";
    exit(1);
  }
  
  if(loc_muntz_seq.size()==1)
  {
    if(loc_muntz_seq[0]>0)
    {
      lambda_min = static_cast<float128>(0.0);
      lambda_max = lambda_min = static_cast<float128>(loc_muntz_seq[0]);
    }
    else
    {
      lambda_min = static_cast<float128>(loc_muntz_seq[0]);
      lambda_max = static_cast<float128>(0.0);
    }
  }
  else
  {
    lambda_min = static_cast<float128>(loc_muntz_seq[0]);
    lambda_max = static_cast<float128>(loc_muntz_seq.back());
  }
  
  if(loud_mode)
  {
    // Print on-screen the inputs recap information
    std::cout << "\n ** Accepted sequence of exponents ** \n";
    std::cout << std::setprecision(std::numeric_limits<float>::max_digits10)
              << "    {" << muntz_sequence[0];
    for(int k=1; k < muntz_sequence.size(); k++)
    {
      std::cout << std::setprecision(std::numeric_limits<float>::max_digits10)
                << ", " << muntz_sequence[k];
    }
    std::cout << "}";
    std::cout << std::setprecision(std::numeric_limits<float>::max_digits10)
              << "\n ** Lambda_min = " << lambda_min 
              << ", Lambda_max = " << lambda_max << " **"
              << std::endl;
  }
    
  // Compute, or stream through, the minimum number of quadrature nodes n_min
  if(compute_n_min)
  {
    n_min = computeNumNodes(lambda_min, lambda_max);
  }
  else
  {
    n_min = num_nodes;
  }

  // Generate and return the outputs
  std::vector<float128> lambdas = {lambda_min, lambda_max};
  return std::make_tuple(n_min, lambdas);
}
template std::tuple<int, std::vector<float128>> manageData<float128>(std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence);
template std::tuple<int, std::vector<float128>> manageData<double>(std::vector<double>& muntz_sequence, std::vector<double>& coeff_sequence);

/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: [n_min, {beta_min, beta_max}] = streamMonMapData(n)
//                
//          INPUT: - n = output of function 'manageData'
//
//         OUTPUT: - n_min = minimum possible number of nodes listed in the 
//                           'data/TabulatedErrorValues.csv' file to obtain 
//                           double-precise numerical integration
//                 - beta_min = minimum value for the exponent of the post-map polynomial
//                 - beta_max = maximum value for the exponent of the post-map polynomial
//
//    DESCRIPTION: the monomial transformation gamma: (0,1) -> (0,1) is uniquely 
//                 identified by its order r which in turn requires the knowledge of
//                 beta_min/beta_max. This method scans the tabulated vales in the 
//                 'data/TabulatedErrorValues.csv' file to extract such values
//                 according to the input n (number of quadrature samples).
//
/////////////////////////////////////////////////////////////////////////////////////////

std::tuple<int, std::vector<float128>, int> streamMonMapData(const int& comp_num_nodes)
{
  int n_min;
  float128 beta_min, beta_max;
  extern bool loud_mode;

  int it = 0;
  while(true)
  {
    std::vector<std::string> _betas = BETAS[it];
    n_min = stoi(_betas[0]);

    if(n_min >= comp_num_nodes)
    {
      std::string _beta_min = _betas[1];
      beta_min = static_cast<float128>(_beta_min);
      std::string _beta_max = _betas[2];
      beta_max = static_cast<float128>(_beta_max);
      if(loud_mode)
      {
        std::cout << std::setprecision(std::numeric_limits<float>::max_digits10)
                  << " ――――――――――――――――――――――――――――――――――――――――――――――――――"
                  << "\n ** N_min = " << n_min
                  << "\n ** Beta_min = " << beta_min
                  << ", Beta_max = " << beta_max
                  << " **";
      }

      break;
    }
    else {
        it++;
    }
  }
  std::vector<float128> betas = {beta_min, beta_max};
  return std::make_tuple(n_min, betas, it);
}


/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: optimiseData(quadrature_parameters, muntz_sequence, coeff_sequence)
//                
//          INPUT: - quadrature_parameters = output of function 'computeQuadParams'
//                 - muntz_sequence = sequence of real exponents of the polynomial
//                 - coeff_sequence = sequence of real coefficients of the polynomial
//
//         OUTPUT: no outputs
//
//    DESCRIPTION: the method automatically selects the most optimised format possible
//                 (between double and float128) to export the output results, i.e. the
//                 transformed (new) quadrature nodes and weights, to reach the 
//                 prescribed relative precision for the integration (e.g. double 
//                 precision). The optimality here is meant as the f.p. format with lowest
//                 precision that still allows to retain a machine-epsilon accuracy (we 
//                 specified double-precision however the procedure is easily generalised
//                 using the method's templatisation) for the relative error of the
//                 integral computed using the monomial transformation quadrature rule.
//
/////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
void optimiseData(std::tuple<std::vector<float128>, std::vector<float128>, std::vector<float128>, std::vector<float128>>& quad_params, std::vector<T>& muntz_sequence, std::vector<T>& coeff_sequence)
{
  const double ACC = 2*EPS;
  bool print_primitive = false;

  extern bool loud_mode;
  
  float128 I34_old, I34_new = computeQuadRule(std::get<0>(quad_params), std::get<1>(quad_params), muntz_sequence, coeff_sequence);
  float128 E34_old, E34_new = computeExactError(I34_new, muntz_sequence, coeff_sequence, print_primitive);

  // Degrade quadrature parameters to double
  std::vector<double> double_nodes = castVector(std::get<0>(quad_params), std::numeric_limits<double>::epsilon());
  std::vector<double> double_weights = castVector(std::get<1>(quad_params), std::numeric_limits<double>::epsilon());
  // Compute I_n(f) and E_n(f) with double precise quadrature parameters
  float128 I16_old, I16_new = computeQuadRule(double_nodes, double_weights, muntz_sequence, coeff_sequence);
  float128 E16_old, E16_new = computeExactError(I16_new, muntz_sequence, coeff_sequence, print_primitive);

  if(fabs(E34_new - E16_new) <= ACC)
  {// Nodes and weights succefully optimised float128 -> double
    if(loud_mode)
    {
      std::cout << " ――――――――――――――――――――――――――――――――――――――――――――――――――"
                << "\n ** Using double f.p. format for nodes and weights **"
                << std::endl;
      print_primitive = true;
      // Degrade classical G-L parameters with double format
      std::vector<double> old_nodes = castVector(std::get<2>(quad_params), std::numeric_limits<double>::epsilon());
      std::vector<double> old_weights = castVector(std::get<3>(quad_params), std::numeric_limits<double>::epsilon());
      // Compute classical G-L quadrature with double format parameters
      I16_old = computeQuadRule(old_nodes, old_weights, muntz_sequence, coeff_sequence);
      E16_old = computeExactError(I16_old, muntz_sequence, coeff_sequence, print_primitive);
      // Print on-screen quadrature results
      std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
                << " **\n ** I_n(p(x)) = "
                << I34_new
                << "  [with parameters in float128 precision]"
                << " **\n ** E_n(p(x)) = "
                << E34_new
                << "  [with parameters in float128 precision]"
                << " **\n ** I_n(p(x)) = "
                << I16_new
                << "  [with parameters in double precision]"
                << " **\n ** E_n(p(x)) = " 
                << E16_new
                << " [with parameters in double precision] **" << std::endl;
      // Check on asbolute accuracy of the relative error
      if(E16_new >= 1e1*ACC)
      {
        std::cout << "\n ** WARNING: the input polynomial is integrated with less than double-precision **" << std::endl;
      }
    }
    // Export all the results {In, En, I, E} in double
    exportNewData(double_nodes, double_weights, {I16_new, E16_new, I16_old, E16_old});
    // Print closing message and return
    if(loud_mode){std::cout << "\n\n ** MTQR HAS TERMINATED **\n";}
    return;
  }
  if(loud_mode)
  {
    std::cout << " ――――――――――――――――――――――――――――――――――――――――――――――――――"
              << "\n ** Using quadruple f.p. format for nodes and weights **"
              << std::endl;
    print_primitive = true;
    // Compute classical G-L quadrature with float128 format parameters
    I34_old = computeQuadRule(std::get<2>(quad_params), std::get<3>(quad_params), muntz_sequence, coeff_sequence);
    E34_old = computeExactError(I34_old, muntz_sequence, coeff_sequence, print_primitive);
    // Print on-screen quadrature results
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
              << " **\n ** I_n(p(x)) = "
              << I34_new
              << "  [with parameters in float128 precision]"
              << " **\n ** E_n(p(x)) = "
              << E34_new
              << "  [with parameters in float128 precision]"
              << " **\n ** I_n(p(x)) = "
              << I16_new
              << "  [with parameters in double precision]"
              << " **\n ** E_n(p(x)) = " 
              << E16_new
              << " [with parameters in double precision] **" << std::endl;
    // Check on asbolute accuracy of the relative error
    if(E16_new <= 1e2*ACC)
    {
      std::cout << "\n ** WARNING: the input polynomial is integrated with less than double-precision **" << std::endl;
    }
  }
  // Export all the results {In, En, I, E} in float128
  exportNewData(std::get<0>(quad_params), std::get<1>(quad_params), {I34_new, E34_new, I34_old, E34_old});
  // Print closing message
  if(loud_mode){std::cout << "\n\n ** MTQR HAS TERMINATED **\n";}
}
template void optimiseData(std::tuple<std::vector<float128>, std::vector<float128>, std::vector<float128>, std::vector<float128>>& quad_params, std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence); 
template void optimiseData(std::tuple<std::vector<float128>, std::vector<float128>, std::vector<float128>, std::vector<float128>>& quad_params, std::vector<double>& muntz_sequence, std::vector<double>& coeff_sequence); 


/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: exportNewData(optim_nodes, optim_weights, [I_new, E_new, I_old, E_old])
//                
//          INPUT: - optim_nodes = output of function 'computeQuadParams' optimised by
//                                 'optimiseData'
//                 - optim_weights = output of function 'computeQuadParams' optimised by
//                                 'optimiseData'
//                 - [I_new, E_new, I_old, E_old] = values of the numerical quadratures
//                                  (I_new, I_old) and relative errors (E_new, E_old)
//                                  obtained using the nodes and weights provided 
//                                  by the monomial transformation and classical G-L
//                                  quadrature rules respectively.
//
//         OUTPUT: no outputs
//
//    DESCRIPTION: once the monomial transformation quadrature rule is completed, data is
//                 streamed in output and written to the appropriate text files. The data
//                 is organised in three separate files:
//                    - 'Results.txt' containins recap informations about the monomial
//                      transformation quadrature rule such as the parameters of the map
//                      (gamma), number of quadrature samples (n_min), etc...
//                    - 'Nodes.txt' listing the new nodes computed by the quadrature rule
//                    - 'Weights.txt' listing the new weights of the above rule
//
/////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
void exportNewData(const std::vector<T>& nodes, const std::vector<T>& weights, const std::vector<float128>& output_data)
{
  extern bool loud_mode;
  if(loud_mode)
  {
    // Write Results.txt file containing transformation review and integral outputs
    std::ofstream Results;
    Results.open("output/Results.txt", std::ios_base::app);
    Results << std::setprecision(std::numeric_limits<float128>::max_digits10)
            << "\n        I_n(f) = "
            << output_data[0]
            << ", E_n(f) = "
            << output_data[1]
            << "  (monomial transformation quadrature rule)"
            << "\n        I_n(f) = "
            << output_data[2]
            << ", E_n(f) = "
            << output_data[3]
            << "  (classical G-L quadrature rule)";
    Results.close();

    // Write Nodes.txt file containing the monomial transformation quadrature rule's nodes
    std::ofstream Nodes;
    Nodes.open("output/Nodes.txt");
    for(int k = 0; k <= nodes.size() - 1; k++)
    {
      Nodes << std::setprecision(std::numeric_limits<T>::max_digits10)
            << nodes[k] << std::endl;
    }
    Nodes.close();

    // Write Weights.txt file containing the monomial transformation quadrature rule's weights
    std::ofstream Weights;
    Weights.open("output/Weights.txt");
    for(int k = 0; k <= weights.size() - 1; k++)
    {
      Weights << std::setprecision(std::numeric_limits<T>::max_digits10)
              << weights[k] << std::endl;
    }
    Weights.close();
  }
  else
  {
    std::ofstream Outputs;
    Outputs.open("outputs.txt", std::ios_base::app);
    Outputs << std::setprecision(std::numeric_limits<float128>::max_digits10)
            << "\n    E_n(f) = "
            << output_data[1]
            << "\n\nNodes:" << std::endl;
    for(int k = 0; k <= nodes.size() - 1; k++)
    {
      Outputs << std::setprecision(std::numeric_limits<T>::max_digits10)
              << "    "
              << nodes[k] << std::endl;
    }
    Outputs << "\nWeights:" << std::endl;
    for(int k = 0; k <= nodes.size() - 1; k++)
    {
      Outputs << std::setprecision(std::numeric_limits<T>::max_digits10)
              << "    "
              << weights[k] << std::endl;
    }
  }
}
template void exportNewData<float128>(const std::vector<float128>& nodes, const std::vector<float128>& weights, const std::vector<float128>& output_data); // Template mock instantiation for non-inline function
template void exportNewData<double>(const std::vector<double>& nodes, const std::vector<double>& weights, const std::vector<float128>& output_data); // Template mock instantiation for non-inline function