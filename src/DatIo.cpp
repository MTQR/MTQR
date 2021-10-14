//---------------------------------------------------------------------------------------
// File:      src/DatIo.cpp
//
// Library:   QUASIMONT-QUAdrature of SIngular polynomials using MONomial Transformations:
//                      a C++ library for high precision integration of singular 
//                      polynomials of non-integer degree
//
// Authors:   Guido Lombardi, PhD, Davide Papapicco
//
// Institute: Politecnico di Torino
//            C.so Duca degli Abruzzi, 24 - Torino (TO), Italia
//            Department of Electronics and Telecommunications (DET)
//            Electromagnetic modelling and applications Research Group
//---------------------------------------------------------------------------------------

#include "Quasimont.h"

/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: generateTabData()
//                
//          INPUT: no inputs
//
//         OUTPUT: no ouptus
//
//    DESCRIPTION: the entirety of the library relies on tabulated values of beta_min 
//                 and beta_max for each even value of number of (quadrature) nodes 
//                 and for this reason it is shipped with those values in  the
//                 'data/TabulatedErrorValues.csv' file. However should the file be 
//                 corrupted or deleted this routine, when triggere by the absence of 
//                 the file itself, reconstructs such table from scratch.
//
/////////////////////////////////////////////////////////////////////////////////////////

void generateTabData()
{

  std::ofstream data("../data/TabulatedErrorValues.csv");

  int num_samples = 15000, prev_j;
  float128 start = -1, end = 1300;
  std::vector<float128> lambda = linspace(start, end, num_samples);

  float128 tab_beta_min, tab_beta_max;
  float1k exact;

  double dp_epsilon = std::numeric_limits<double>::epsilon();

  std::string envelope = "yes";

  std::cout << "\nTabulated values of beta_min/max not found!\nComputing as function of n (number of nodes)\n";

  int min_nodes = 5, max_nodes = 50;

  for(int n = min_nodes; n <= max_nodes; n++)
  {

    int nodes = 2*n;
    std::cout << "\nComputing for n = " << nodes << std::endl;

    for(int k = 0; k < num_samples - 2; k++)
    {

      std::cout << "    beta_min = " << lambda[k+1] << "\r";

      exact = computeEstimate(lambda[k+1], nodes, envelope);

      if(exact <= dp_epsilon)
      {
        std::cout << std::endl;
        tab_beta_min = lambda[k+1];

        if(prev_j < k+1)
        {
          prev_j = k+1;
        }

        for(int j = prev_j+1; j < num_samples - 2; j++)
        {
          std::cout << "    beta_max = " << lambda[j] << "\r";

          exact = computeEstimate(lambda[j], nodes, envelope);

          if(exact >= dp_epsilon)
          {
            tab_beta_max = lambda[j];
            prev_j = j;

            break;
          }
        }

        break;
      }
    }

    data << std::setprecision(std::numeric_limits<float1k>::max_digits10) << nodes << "," << tab_beta_min << "," << tab_beta_max << "\n";
  }

  data.close();

  std::cout << "\n\nComputation completed!\nTabulated data is now available in the 'HIPER-SINK/data' subdirectory.\nProgram can now be executed!\n" << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: [n, {lambda_min, lambda_max}] = getInputData(muntz_sequence, 
//                                                              coeff_sequence)
//                
//        INPUT: - muntz_sequence = sequence of real exponents of the polynomial
//               - coeff_sequence = sequence of real coefficients of the polynomial
//
//       OUTPUT: - n = output of function 'computeNumNodes'
//               - lambda_min = minimum exponent in the input "muntz_sequence"
//               - lambda_max = maximum exponent in the input "muntz_sequence"
//
//    DESCRIPTION: the user-input polynomial is provided to the library via Main.cpp; the
//                 polynomial itself is specified via a unordered sequence of 
//                 coefficients and exponents of the various monomials in the polynomial.
//                 Once those input are read, checks have to be made in order to validate
//                 the proper functioning of the library; those are:
//                    - the number of exponents and the number of coefficients coincide;
//                    - the input polynomial is at least a binomial (otherwise the 
//                      routine further CLI user-input is required, see lines 86~93 in
//                      the 'src/MonMap.cpp' file); 
//                    - lambda_min > -1 (otherwise the program exits);
//                 Once those checks are ran the exponents' sequence is sorted locally
//                 and lambda_min/lambda_max are thus identified and outputted alongside
//                 the associated number of nodes computed by the function 
//                 'computeNumNodes' (see lines 110~175 in the 'src/MonMap.cpp' file).
//
/////////////////////////////////////////////////////////////////////////////////////////

std::tuple<int, std::vector<float128>> getInputData(std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence)
{

  std::cout << std::endl;
  std::cout << "    |――――――――――――――――――――――――――――――――――|\n"
        << "    |  ** MONOMIAL QUADRATURE RULE **  |\n"
        << "    |――――――――――――――――――――――――――――――――――|\n";

  if(muntz_sequence.size()==coeff_sequence.size())
  {
    std::cout << "\n\n Input polynomial p(x) = ";
    for(int k=0; k < muntz_sequence.size(); k++)
    {
      if(coeff_sequence[k]>0)
      {
        std::cout << "+";
      }

      if(coeff_sequence[k]==1)
      {
        std::cout << "x^(" << muntz_sequence[k] << ") ";
      }
      else
      {
        std::cout << coeff_sequence[k] << "*x^(" << muntz_sequence[k] << ")";
      }
    }

    std::cout << std::endl;
  }
  else
  {
    std::cout << "\n   ** ERROR ** The number of exponents doesn't match the number of coefficients in the input polynomial\n";
  }

  int num_nodes, n_min;
  float128 lambda_max, lambda_min;
  bool compute_n_min = true;

  if(muntz_sequence.size()==1)
    {
      if(muntz_sequence[0] <= -1.0)
      {
        std::cout << "\n   ** ERROR ** Lambda_min has to be (strictly) greater than -1 to be in a Muntz sequence\n";
      }

      float128 additional_lambda;

      std::cout << "\n   ** WARNING ** Your input is a monomial of non-integer degree"
            << "\n                 You need at least a binomial to achieve double-precision quadrature"
                << "\n                 How do you want to proceed?"
                << "\n   [enter 'nodes' to specify the number of nodes or 'lambda' for the exponent]\n"
            << "\n                 Input: ";

      std::string input;

      std::cin >> input;

      if(input.compare("nodes") == 0)
      {
        std::cout << "\n\nPlease specify the desired number of quadrature nodes (number must be even): ";
        std::cin >> num_nodes;

        compute_n_min = false;

        additional_lambda = computeLambda(muntz_sequence[0], num_nodes); // lambda_min = muntz_sequence[0]
      }
      else
      {
        std::cout << "\n\nPlease enter the exponent value with '.' separating the decimal digits from the integer part [the more decimal digits the better the precision is]: ";
        std::cin >> additional_lambda;
      }

      muntz_sequence.push_back(additional_lambda);
      coeff_sequence.push_back(1.0);
    }

  std::vector<float128> loc_muntz_seq = muntz_sequence;

  sort(loc_muntz_seq.begin(), loc_muntz_seq.end());

  if(loc_muntz_seq[0] <= -1.0)
  {
    std::cout << "\n   ** ERROR ** Lambda_min has to be (strictly) greater than -1 to be in a Muntz sequence\n";
  }
  
  if(loc_muntz_seq.size()==1)
  {
    if(loc_muntz_seq[0]>0)
    {
      lambda_min = 0;
      lambda_max = lambda_min = loc_muntz_seq[0];
    }
    else
    {
      lambda_min = loc_muntz_seq[0];
      lambda_max = 0;
    }
  }
  else
  {
    lambda_min = loc_muntz_seq[0];
    lambda_max = loc_muntz_seq.back();
  }
  
  std::cout << "\n ** Accepted sequence of exponents ** \n";
  std::cout << "    {" << muntz_sequence[0];
  for(int k=1; k < muntz_sequence.size(); k++)
  {
    std::cout << ", " << muntz_sequence[k];
  }
  std::cout << "}";

  std::cout << "\n ** Lambda_min = " << lambda_min << ", Lambda_max = " << lambda_max << " **" << std::endl;

  std::vector<float128> lambdas = {lambda_min, lambda_max};

  if(compute_n_min)
  {
    n_min = computeNumNodes(lambda_min, lambda_max);
  }
  else
  {
    n_min = num_nodes;
  }

  return std::make_tuple(n_min, lambdas);
}

/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: [n_min, {beta_min, beta_max}] = retrieveMonData(n)
//                
//        INPUT: - n = output of function 'getInputData'
//
//       OUTPUT: - n_min = minimum possible (even) number of nodes from the 
//                         'data/TabulatedErrorValues.csv' file
//               - beta_min = minimum value for the exponent of the post-map polynomial
//               - beta_max = maximum value for the exponent of the post-map polynomial
//
//    DESCRIPTION: the monomial transformation gamma: [0,1] -> [0,1] is uniquely 
//                 identified by its order r which in turn requires the knowledge of
//                 beta_min/beta_max, alongside lambda_min/lambda_max, to be computed
//                 (see lines 178~211 in the 'src/MonMap.cpp' file). This method scans
//                 the tabulated vales in the 'data/TabulatedErrorValues.csv' file to
//                 extract the beta_min/beta_max and n_min required by the monomial
//                 quadrature rule according to the specified number of nodes as either
//                 computed by the function 'computeNumNodes' (see line 243) or provided
//                 as user-input (see lines 168~199, 247).
//
/////////////////////////////////////////////////////////////////////////////////////////

std::tuple<int, std::vector<double>> retrieveMonData(const int& comp_num_nodes)
{
  int n_min;
  double beta_min, beta_max;

  std::string line, column;

  std::ifstream datafile;
  datafile.open("../data/TabulatedErrorValues.csv");

  while(std::getline(datafile, line))
  {
    
    std::stringstream column_string(line);

    std::vector<std::string> row;

    while(std::getline(column_string, column, ','))
    {
      row.push_back(column);
    }

    n_min = stoi(row[0]);

    if(n_min >= comp_num_nodes)
    {

      beta_min = stof(row[1]);
      beta_max = stof(row[2]);

      std::cout << " ――――――――――――――――――――――――――――――――――――――――――――――――――"
                << "\n ** Beta_min = " << row[1].substr(0,10)
                << ", Beta_max = " << row[2].substr(0,10)
                << " **";

      break;
    }

  }

  std::vector<double> betas = {beta_min, beta_max};

  return std::make_tuple(n_min, betas);
}

/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: degradeData()
//                
//          INPUT: - [J, {new_x, new_w, old_x, old_w}] = output of function
//                                                       'computeParams'
//                 - muntz_sequence = sequence of real exponents of the polynomial
//                 - coeff_sequence = sequence of real coefficients of the polynomial
//
//         OUTPUT: no outputs
//
//    DESCRIPTION: TBA
//
/////////////////////////////////////////////////////////////////////////////////////////


void degradeData(std::tuple<double, std::vector<float50>, std::vector<float50>, std::vector<float50>, std::vector<float50>>& quad_params, std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence, const std::vector<double>& interval)
{
  // jacobian = std::get<0>(quad_params), new_nodes = std::get<1>(quad_params), new_weights = std::get<2>(quad_params)
  float50 float50_quadrature = std::get<0>(quad_params)*computeQuadGl(std::get<1>(quad_params), std::get<2>(quad_params), muntz_sequence, coeff_sequence);
  float50 float50_error = computeError(float50_quadrature, muntz_sequence, coeff_sequence, interval);

  double epsilon = std::numeric_limits<double>::epsilon();

  std::string data_type = "float50";

  if(float50_error <= epsilon)
  {

    std::vector<float128> float128_nodes, float128_weights; 
    castVector(std::get<1>(quad_params), float128_nodes);
    castVector(std::get<2>(quad_params), float128_weights);
      
      float128 float128_quadrature = std::get<0>(quad_params)*computeQuadGl(float128_nodes, float128_weights, muntz_sequence, coeff_sequence);
      float128 float128_error = computeError(float128_quadrature, muntz_sequence, coeff_sequence, interval);

      if(float128_error <= epsilon && float128_error != 0)
      {

        data_type = "float128";

        std::vector<double> double_nodes, double_weights;
        castVector(std::get<1>(quad_params), double_nodes); 
        castVector(std::get<2>(quad_params), double_weights);
        
        double double_quadrature = static_cast<double>(std::get<0>(quad_params)*computeQuadGl(double_nodes, double_weights, muntz_sequence, coeff_sequence));
        double double_error = static_cast<double>(computeError(double_quadrature, muntz_sequence, coeff_sequence, interval));

        if(double_error <= epsilon && double_error != 0)
        {

          data_type = "double";

          std::vector<double> old_nodes, old_weights;
          castVector(std::get<3>(quad_params), old_nodes);
          castVector(std::get<4>(quad_params), old_weights);
          
          double classical_quadrature = static_cast<double>(std::get<0>(quad_params)*computeQuadGl(old_nodes, old_weights, muntz_sequence, coeff_sequence));
          double classical_error = static_cast<double>(computeError(classical_quadrature, muntz_sequence, coeff_sequence, interval));
          
          std::vector<double> output_data = {double_quadrature, double_error, classical_quadrature, classical_error};
          exportData(quad_params, output_data, data_type);
        }
        else
        {
          std::vector<float128> old_nodes, old_weights;
          castVector(std::get<3>(quad_params), old_nodes);
          castVector(std::get<4>(quad_params), old_weights);

          float128 classical_quadrature = std::get<0>(quad_params)*computeQuadGl(old_nodes, old_weights, muntz_sequence, coeff_sequence);
          float128 classical_error = computeError(classical_quadrature, muntz_sequence, coeff_sequence, interval);
          
          std::vector<float128> output_data = {float128_quadrature, float128_error, classical_quadrature, classical_error};
          exportData(quad_params, output_data, data_type);
        }
      }
      else
      {
        // jacobian = std::get<0>(quad_params), old_nodes = std::get<3>(quad_params), old_weights = std::get<4>(quad_params)
        float50 classical_quadrature = std::get<0>(quad_params)*computeQuadGl(std::get<3>(quad_params), std::get<4>(quad_params), muntz_sequence, coeff_sequence);
        float50 classical_error = computeError(classical_quadrature, muntz_sequence, coeff_sequence, interval);

        std::vector<float50> output_data = {float50_quadrature, float50_error, classical_quadrature, classical_error};
        exportData(quad_params, output_data, data_type);
      }
  }
  else
  {
    // jacobian = std::get<0>(quad_params), old_nodes = std::get<3>(quad_params), old_weights = std::get<4>(quad_params)
    float50 classical_quadrature = std::get<0>(quad_params)*computeQuadGl(std::get<3>(quad_params), std::get<4>(quad_params), muntz_sequence, coeff_sequence);
    float50 classical_error = computeError(classical_quadrature, muntz_sequence, coeff_sequence, interval);

    std::vector<float50> output_data = {float50_quadrature, float50_error, classical_quadrature, classical_error};
    exportData(quad_params, output_data, data_type);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: exportData({lambda_min, lambda_max}, [n_min, {beta_min, beta_max}],
//                            [], {post_map_integral, pre_map_integral}, r)
//                
//          INPUT: - {lambda_min, lambda_max} = output of function 'getInputData'
//                 - [n_min, {beta_min, beta_max}] = output of function 'retrieveMonData'
//                 - [] = output of function 'computeParams'
//                 - {post_map_integral, pre_map_integral} = output of function
//                                                           'computeQuadGl'
//                 - r = output of function 'computeOrder'
//
//         OUTPUT: no outputs
//
//    DESCRIPTION: once the monomial quadrature rule is completed the resulting data is
//                 streamed in output text files inside the output subdirectory created,
//                 by this routine, within the calling directory of the executable of 
//                 this library. The output data is splitted in three files with
//                 'Results.txt' containing recap informations of the execution 
//                 (including the resulting integral) and 'Nodes.txt' and 'Weights.txt'
//                 containg the classic and new G-L nodes and weights respectively.
//
/////////////////////////////////////////////////////////////////////////////////////////


template<typename type>
void exportData(const std::tuple<double, std::vector<float50>, std::vector<float50>, std::vector<float50>, std::vector<float50>>& quad_params, const std::vector<type>& output_data, const std::string& data_type)
{
  std::cout << " ――――――――――――――――――――――――――――――――――――――――――――――――――"
        << "\n ** Using " << data_type << " format for nodes and weights"
        << std::setprecision(std::numeric_limits<float128>::max_digits10)
        << "\n ** I_n(p(x)) = "
        << output_data[0]
        << " **\n ** E_n(p(x)) = " 
        << output_data[1]
        << " **" << std::endl;

  // jacobian = std::get<0>(quad_params), classical_nodes = std::get<3>(quad_params), classical_weights = std::get<4>(quad_params);

  std::vector<type> new_nodes, new_weights;
  castVector(std::get<1>(quad_params), new_nodes);
  castVector(std::get<2>(quad_params), new_weights);
  
  std::vector<type> old_nodes, old_weights;
  castVector(std::get<3>(quad_params), old_nodes);
  castVector(std::get<4>(quad_params), old_weights);

  std::ofstream results;
  results.open("output/Results.txt", std::ios_base::app);

  results << std::setprecision(std::numeric_limits<float128>::max_digits10)
          << "\n        I_n(f) = "
          << output_data[0]
          << ", E_n(f) = "
          << output_data[1]
          << "  (monomial quadrature rule)"
          << "\n        I_n(f) = "
          << output_data[2]
          << ", E_n(f) = "
          << output_data[3]
          << "  (classical G-L rule)";

    results.close();

    // Export G-L nodes (old and new)

    std::string nodes_file = "output/Nodes.txt";

  std::ofstream Nodes;
  Nodes.open(nodes_file);

    Nodes << "\nTRANSFORMED G-L NODES (MONOMIAL MAPPED):"; 

  for(int k = 0; k <= new_nodes.size() - 1; k++)
  {

    Nodes << std::setprecision(std::numeric_limits<float50>::max_digits10)
        << "\n        x_"
        << k
        << ": "
        << new_nodes[k];

  }

  Nodes << "\n\nCLASSIC G-L NODES (AFFINE MAPPED):"; 

  for(int k = 0; k <= old_nodes.size() - 1; k++)
  {

    Nodes << std::setprecision(std::numeric_limits<float50>::max_digits10)
        << "\n        x_"
        << k
        << ": "
        << old_nodes[k];

  }

  Nodes.close();

  // Export G-L weights (old and new)

    std::string weights_file = "output/Weights.txt";

  std::ofstream Weights;
  Weights.open(weights_file);

    Weights << "\nTRANSFORMED G-L WEIGHTS (MONOMIAL MAPPED):"; 

  for(int k = 0; k <= new_weights.size() - 1; k++)
  {

    Weights << std::setprecision(std::numeric_limits<float50>::max_digits10)
        << "\n        w_"
        << k
        << ": "
        << new_weights[k];

  }

  Weights << "\n\nCLASSIC G-L WEIGHTS (AFFINE MAPPED):"; 

  for(int k = 0; k <= old_weights.size() - 1; k++)
  {

    Weights << std::setprecision(std::numeric_limits<float50>::max_digits10)
        << "\n        w_"
        << k
        << ": "
        << old_weights[k];

  }

  Weights.close();
}
template void exportData<float50>(const std::tuple<double, std::vector<float50>, std::vector<float50>, std::vector<float50>, std::vector<float50>>& quad_params, const std::vector<float50>& output_data, const std::string& data_type); // Template mock instantiation for non-inline function
template void exportData<float128>(const std::tuple<double, std::vector<float50>, std::vector<float50>, std::vector<float50>, std::vector<float50>>& quad_params, const std::vector<float128>& output_data, const std::string& data_type); // Template mock instantiation for non-inline function
template void exportData<double>(const std::tuple<double, std::vector<float50>, std::vector<float50>, std::vector<float50>, std::vector<float50>>& quad_params, const std::vector<double>& output_data, const std::string& data_type); // Template mock instantiation for non-inline function