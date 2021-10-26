//---------------------------------------------------------------------------------------
// File:      src/DatIo.cpp
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

void checkTabData()
{
  std::ifstream datafile;
  datafile.open("../data/TabulatedErrorValues.csv");

  if(datafile)
  {}// Tabulated data exists ==> Pass 
  else
  {
    std::ofstream data("../data/TabulatedErrorValues.csv");

    int num_samples = 15000, prev_j;
    float128 start = -1, end = 1300;
    std::vector<float128> lambda = linspace(start, end, num_samples);

    float128 tab_beta_min, tab_beta_max;
    float1k exact;

    std::string envelope = "yes";

    std::cout << "\n ** Tabulated values of beta_min/max not found ** \n ** Computing as function of n (number of nodes) ** \n";

    int min_nodes = 5, max_nodes = 50;

    for(int n = min_nodes; n <= max_nodes; n++)
    {
      int nodes = 2*n;
      std::cout << "\nComputing for n = " << nodes << std::endl;

      for(int k = 0; k < num_samples - 2; k++)
      {
        std::cout << "    beta_min = " << lambda[k+1] << "\r";

        exact = computeErrorEstimate(lambda[k+1], nodes, envelope);

        if(exact <= EPS)
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

            exact = computeErrorEstimate(lambda[j], nodes, envelope);

            if(exact >= EPS)
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

    std::cout << "\n\n ** Tabulated data is now available in 'QUASIMONT/data' subdirectory **\n ** Proceeding with the application **\n" << std::endl;;
  }
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

template <typename type>
std::tuple<int, std::vector<float128>> manageData(std::vector<type>& muntz_sequence, std::vector<type>& coeff_sequence)
{
  std::cout << std::endl;
  std::cout << "    |――――――――――――――――――――――――――――――――――|\n"
            << "    |          ** QUASIMONT **         |\n"
            << "    |  ** MONOMIAL QUADRATURE RULE **  |\n"
            << "    |――――――――――――――――――――――――――――――――――|\n";

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

  int num_nodes, n_min;
  float128 lambda_max, lambda_min;
  bool compute_n_min = true;

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
                << "\n               QUASIMONT needs a binomial for double-precision quadrature."
                << "\n               How do you proceed? ['nodes' for n_min ~ 'lambda' for lambda_max]"
                << "\n               Input: ";
      std::cin >> input;

      if(input.compare("nodes") == 0)
      {
        std::cout << "\n\nPlease specify the desired number of quadrature nodes (number must be even): ";
        std::cin >> num_nodes;

        compute_n_min = false;

        float128 local_lambda_min = static_cast<float128>(muntz_sequence[0]);
        additional_lambda = computeLambdaMax(local_lambda_min, num_nodes);
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

      muntz_sequence.push_back(static_cast<type>(additional_lambda));
      coeff_sequence.push_back(1.0);
    }

  // Sorting the input Muntz sequance to extract lambda_min and lambda_max
  std::vector<type> loc_muntz_seq = muntz_sequence;
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
      lambda_min = 0;
      lambda_max = lambda_min = static_cast<float128>(loc_muntz_seq[0]);
    }
    else
    {
      lambda_min = static_cast<float128>(loc_muntz_seq[0]);
      lambda_max = 0;
    }
  }
  else
  {
    lambda_min = static_cast<float128>(loc_muntz_seq[0]);
    lambda_max = static_cast<float128>(loc_muntz_seq.back());
  }
  
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
  
  // Compute or stream through the minimum number of quadrature nodes n_min
  if(compute_n_min)
  {
    n_min = computeNumNodes(lambda_min, lambda_max);
  }
  else
  {
    n_min = num_nodes;
  }
  return std::make_tuple(n_min, {lambda_min, lambda_max});
}
template std::tuple<int, std::vector<float128>> manageData<float50>(std::vector<float50>& muntz_sequence, std::vector<float50>& coeff_sequence);
template std::tuple<int, std::vector<float128>> manageData<float128>(std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence);
template std::tuple<int, std::vector<float128>> manageData<double>(std::vector<double>& muntz_sequence, std::vector<double>& coeff_sequence);


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

std::tuple<int, std::vector<double>> streamMonMapData(const int& comp_num_nodes)
{
  int n_min;
  double beta_min, beta_max;

  std::ifstream datafile;
  datafile.open("../data/TabulatedErrorValues.csv");

  std::string line, column;
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
                << "\n ** N_min = " << comp_num_nodes
                << "\n ** Beta_min = " << row[1].substr(0,10)
                << ", Beta_max = " << row[2].substr(0,10)
                << " **";

      break;
    }
  }

  return std::make_tuple(n_min, {beta_min, beta_max});
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

template <typename type>
void optimiseData(std::tuple<std::vector<float50>, std::vector<float50>, std::vector<float50>, std::vector<float50>>& quad_params, std::vector<type>& muntz_sequence, std::vector<type>& coeff_sequence)
{
  const double ACC = 7*EPS;
  std::cout << "\nAccuracy = " << ACC << std::endl;
  // new_nodes = std::get<0>(quad_params)       // old_nodes = std::get<2>(quad_params)
  // new_weights = std::get<1>(quad_params)     // old_weights = std::get<3>(quad_params)
  std::vector<float50> float50_exp, float50_coef;
  castVector(muntz_sequence, float50_exp);
  castVector(coeff_sequence, float50_coef);
  float50 In_50 = computeQuadGl(std::get<0>(quad_params), std::get<1>(quad_params), float50_exp, float50_coef);
  float50 En_50 = computeExactError(In_50, float50_exp, float50_coef);

  std::cout << "\n\nI_n(f) with float50 = "
            << std::setprecision(std::numeric_limits<float50>::max_digits10)
            << In_50
            << "\nE_n(f) with float50 = "
            << En_50 << std::endl;

  if(En_50 <= ACC)
  {
    // Degrade G-L parameters to float128
    std::vector<float128> float128_nodes, float128_weights, float128_exp, float128_coef; 
    castVector(std::get<0>(quad_params), float128_nodes);
    castVector(std::get<1>(quad_params), float128_weights);
    castVector(muntz_sequence, float128_exp);
    castVector(coeff_sequence, float128_coef);
    // Compute I_n(f) and E_n(f) with float128 G-L parameters
    float128 In_34 = computeQuadGl(float128_nodes, float128_weights, float128_exp, float128_coef);
    float128 En_34 = computeExactError(In_34, float128_exp, float128_coef);

    std::cout << "\nI_n(f) with float128 = "
              << std::setprecision(std::numeric_limits<float50>::max_digits10)
              << In_34
              << "\nE_n(f) with float128 = "
              << En_34 << std::endl;

    if(En_34 <= ACC && En_34 != 0)
    {// Nodes and weights succefully optimised float50 -> float128
      
      // Degrade G-L parameters to double
      std::vector<double> double_nodes, double_weights, double_exp, double_coef;
      castVector(std::get<0>(quad_params), double_nodes); 
      castVector(std::get<1>(quad_params), double_weights);
      castVector(muntz_sequence, double_exp);
      castVector(coeff_sequence, double_coef);
      // Compute I_n(f) and E_n(f) with double G-L parameters
      double In_16 = static_cast<double>(computeQuadGl(double_nodes, double_weights, double_exp, double_coef));
      double En_16 = static_cast<double>(computeExactError(In_16, double_exp, double_coef));

      std::cout << "\n\nI_n(f) with double = "
                << std::setprecision(std::numeric_limits<float50>::max_digits10)
                << In_16
                << "\nE_n(f) with double = "
                << En_16 << std::endl;

      if(En_16 <= ACC && En_16 != 0)
      {// Nodes and weights succesfully optimised float128 -> double
        
        // Degrade classical G-L parameters with double format
        std::vector<double> old_nodes, old_weights;
        castVector(std::get<2>(quad_params), old_nodes);
        castVector(std::get<3>(quad_params), old_weights);
        // Compute classical G-L quadrature with double format parameters
        double I_16 = static_cast<double>(computeQuadGl(old_nodes, old_weights, double_exp, double_coef));
        double E_16 = static_cast<double>(computeExactError(I_16, double_exp, double_coef));
        // Print and export all the results {In, En, I, E} in double
        std::cout << " ――――――――――――――――――――――――――――――――――――――――――――――――――"
                  << "\n ** Using double f.p. format for nodes and weights **";
        exportNewData(double_nodes, double_weights, {In_16, En_16, I_16, E_16});
      }
      else
      {
        // Degrade classical G-L parameters with float128 format
        std::vector<float128> old_nodes, old_weights;
        castVector(std::get<2>(quad_params), old_nodes);
        castVector(std::get<3>(quad_params), old_weights);
        // Compute classical G-L quadrature with float128 format parameters
        float128 I_34 = computeQuadGl(old_nodes, old_weights, float128_exp, float128_coef);
        float128 E_34 = computeExactError(I_34, float128_exp, float128_coef);
        // Print and export all the results {In, En, I, E} in float128
        std::cout << " ――――――――――――――――――――――――――――――――――――――――――――――――――"
                  << "\n ** Using float128 f.p. format for nodes and weights **";
        exportNewData(float128_nodes, float128_weights, {In_34, En_34, I_34, E_34});
      }
    }
  }
  else
  {
    // Compute classical G-L quadrature with float128 format parameters
    float50 I_50 = computeQuadGl(std::get<2>(quad_params), std::get<3>(quad_params), float50_exp, float50_coef);
    float50 E_50 = computeExactError(I_50, float50_exp, float50_coef);
    // Print and export all the results {In, En, I, E} in float50
    std::cout << " ――――――――――――――――――――――――――――――――――――――――――――――――――"
              << "\n ** Using float50 f.p. format for nodes and weights **";
    exportNewData(std::get<0>(quad_params), std::get<1>(quad_params), {In_50, En_50, I_50, E_50});
  }

  std::cout << "\n\n ** QUASIMONT HAS TERMINATED **\n";
}
template void optimiseData<float50>(std::tuple<std::vector<float50>, std::vector<float50>, std::vector<float50>, std::vector<float50>>& quad_params, std::vector<float50>& muntz_sequence, std::vector<float50>& coeff_sequence);
template void optimiseData<float128>(std::tuple<std::vector<float50>, std::vector<float50>, std::vector<float50>, std::vector<float50>>& quad_params, std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence);
template void optimiseData<double>(std::tuple<std::vector<float50>, std::vector<float50>, std::vector<float50>, std::vector<float50>>& quad_params, std::vector<double>& muntz_sequence, std::vector<double>& coeff_sequence);


/////////////////////////////////////////////////////////////////////////////////////////
//
//       FUNCTION: exportNewData({lambda_min, lambda_max}, [n_min, {beta_min, beta_max}],
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
void exportNewData(const std::vector<type>& nodes, const std::vector<type>& weights, const std::vector<type>& output_data)
{
  // Print on-screen quadrature results
  std::cout << std::setprecision(std::numeric_limits<float>::max_digits10)
        << "\n ** I_n(p(x)) = "
        << output_data[0]
        << " **\n ** E_n(p(x)) = " 
        << output_data[1]
        << " **" << std::endl;

  // Write Results.txt file containing transformation review and integral outputs
  std::ofstream Results;
  Results.open("output/Results.txt", std::ios_base::app);
  Results << std::setprecision(std::numeric_limits<float128>::max_digits10)
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
  Results.close();

  // Write Nodes.txt file containing transformed G-L nodes using the monomial map
  std::ofstream Nodes;
  Nodes.open("output/Nodes.txt");
  for(int k = 0; k <= nodes.size() - 1; k++)
  {
    Nodes << std::setprecision(std::numeric_limits<type>::max_digits10)
          << nodes[k] << std::endl;
  }
  Nodes.close();

  // Write Weights.txt file containing transformed G-L weights using the monomial map
  std::ofstream Weights;
  Weights.open("output/Weights.txt");
  for(int k = 0; k <= weights.size() - 1; k++)
  {
    Weights << std::setprecision(std::numeric_limits<type>::max_digits10)
            << weights[k] << std::endl;
  }
  Weights.close();
}
template void exportNewData<float50>(const std::vector<float50>& nodes, const std::vector<float50>& weights, const std::vector<float50>& output_data); // Template mock instantiation for non-inline function
template void exportNewData<float128>(const std::vector<float128>& nodes, const std::vector<float128>& weights, const std::vector<float128>& output_data); // Template mock instantiation for non-inline function
template void exportNewData<double>(const std::vector<double>& nodes, const std::vector<double>& weights, const std::vector<double>& output_data); // Template mock instantiation for non-inline function

/*void appendOldData(const std::vector<float50>& GLnodes, const std::vector<float50>& GLweights)
{
  // Write second part of Nodes.txt by appending classical G-L nodes after the affine map in [0,1]
  std::ofstream Nodes;
  Nodes.open("output/Nodes.txt", std::ios_base::app);

  Nodes << "\n\nOLD G-L NODES IN [0,1]:"; 
  for(int k = 0; k <= GLnodes.size() - 1; k++)
  {
    Nodes << std::setprecision(std::numeric_limits<float50>::max_digits10)
        << "\n        x_"
        << k
        << ": "
        << GLnodes[k];
  }
  Nodes.close();

  // Write second part of Weights.txt by appending classical G-L weights after the affine map in [0,1]
  std::ofstream Weights;
  Weights.open("output/Weights.txt", std::ios_base::app);

  Weights << "\n\nOLD G-L WEIGHTS IN [0,1]:"; 
  for(int k = 0; k <= GLweights.size() - 1; k++)
  {
    Weights << std::setprecision(std::numeric_limits<float50>::max_digits10)
        << "\n        w_"
        << k
        << ": "
        << GLweights[k];
  }
  Weights.close();
}*/