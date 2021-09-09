#include "Prototypes.h"

int main(int argc, char** argv)
{

	// CHECK IF TABULATED VALUES EXIST

	std::ifstream datafile;
	datafile.open("../../data/TabulatedErrorValues.csv");

	if(!datafile)
	{
		#include "GenerateTabValues.h"
	}

	if(argc==1) // EXECUTE THE PROGRAM
	{

		// PRINT INITIAL MESSAGE AND TAKES USER INPUTS

		#include "Inputs.h"

		// DERIVE MINIMUN NUMBER OF NODES

		int n_nodes = compute_n(lambda_max, lambda_min);

		// EXTRACT BETA_MIN AND BETA_MAX

		#include "ExtractData.h"

		// COMPUTE THE MONOMIAL TRANSFORMATION ORDER

		double r_max = (1 + beta_max)/(1 + lambda_max), r_min = (1 + beta_min)/(1 + lambda_min);
		double transf_order = (r_min + r_max)/2;
		
	  std::cout << std::left 
	            << std::setw(11) << ""
	            << std::setw(11) << "r = " + std::to_string(transf_order) + ".\n\n";

		// COMPUTE G-L WEIGHTS AND NODES AND EXPORT THEM ALONG WITH OTHER OUTPUTS AND PLOTS

		export_results(lambda_min, lambda_max, n_min, beta_min, beta_max, transf_order, plot, test);

		std::cout << "Output data available in the 'HIPER-SINK/output/n=" << n_min << "' subdirectoy.\n\nProgram terminated!\n" << std::endl;

	}
	else // EXECUTE THE TEST DRIVERS
	{

		std::string test = argv[1];
		
		#include "GetInputData.h"

		int n_nodes = compute_n(lambda_max, lambda_min);

		#include "ExtractData.h"

		double r_max = (1 + beta_max)/(1 + lambda_max), r_min = (1 + beta_min)/(1 + lambda_min);
		double transf_order = (r_min + r_max)/2;

	  std::cout << std::left 
	            << std::setw(11) << ""
	            << std::setw(11) << "r = " + std::to_string(transf_order) + ".\n\n";

		export_results(lambda_min, lambda_max, n_min, beta_min, beta_max, transf_order, plot, test);

		std::cout << "Output data available in the 'HIPER-SINK/output/" << test << "' subdirectoy.\n\nProgram terminated!\n" << std::endl;

	}

	return 0;

}

template<typename type>
int compute_n(const type& lambda_max, const type& lambda_min)
{

	double dp_epsilon = std::numeric_limits<double>::epsilon();

	float1k l_max = static_cast<float1k>(lambda_max), l_min = static_cast<float1k>(lambda_min); 
	constexpr double c0 = -0.0040693, c1 = 0.00041296, d0 = 7.8147, d2 = 0.10123;

	float1k c2 = d2 + d2*l_min;
	float1k c3 = d0 + d0*l_min + l_min - l_max;
	float1k c4 = boost::math::powm1(1+l_max,3) + 1;

	double coeff7 = static_cast<double>(c1*(boost::math::powm1(c2,3)+1));
	double coeff6 = static_cast<double>(c0*(boost::math::powm1(c2,3)+1));
	double coeff5 = static_cast<double>(3*c1*(boost::math::powm1(c2,2) +1)*c3);
	double coeff4 = static_cast<double>(3*c0*(boost::math::powm1(c2,2) + 1)*c3);
	double coeff3 = static_cast<double>(3*c1*c2*(boost::math::powm1(c3,2)+1));
	double coeff2 = static_cast<double>(3*c0*c2*(boost::math::powm1(c3,2)+1));
	double coeff1 = static_cast<double>(c1*(boost::math::powm1(c3,3)+1));
	double coeff0 = static_cast<double>(c1*(boost::math::powm1(c3,3)+1) - c4);

  double coeff[8] = {coeff0, coeff1, coeff2, coeff3, coeff4, coeff5, coeff6, coeff7};
  double n[14];

  gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc	(8);
  gsl_poly_complex_solve(coeff, 8, w, n);
  gsl_poly_complex_workspace_free(w);

  int n_min, flag = 0;

  for (int k = 0; k < 7; k++)
  {
    if(fabs(n[2*k+1]) < dp_epsilon)
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

template<typename type>
float1k compute_estimate(const type& input_lambda, const int& num_nodes, const std::string& envelope)
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

template<typename type>
std::vector<type> linspace(const type& start_type, const type& end_type, const int& num_steps)
{

  std::vector<type> linspaced;

  type start = static_cast<type>(start_type);
  type end = static_cast<type>(end_type);
  type steps = static_cast<type>(num_steps);

  if (steps == 0)
  { 
  	return linspaced;
  }

  if (steps == 1) 
  {
    linspaced.push_back(start);
    return linspaced;
  }

  type delta = (end - start) / (steps - 1);

  for(int k = 0; k < steps - 1; k++)
  {
    linspaced.push_back(start + delta*k);
  }

  linspaced.push_back(end);
  return linspaced;

}

void export_results(const double& lambda_min, const double& lambda_max, const int& n_min, const double& beta_min, const double& beta_max, const double& r, const std::string& plot_flag, const std::string& test_flag)
{

	std::string results_file;

	if(test_flag.compare("n") == 0)
	{
		std::string mkdir = "mkdir -p ";
		std::string output_dir = "../../output/n=";
		std::string results_dir = mkdir + output_dir + std::to_string(n_min);
		results_file = output_dir + std::to_string(n_min) + "/Results.txt";
		system(results_dir.c_str());
	}
	else
	{
		std::string output_dir = "../../output/";
		results_file = output_dir + test_flag + "/Results.txt";
	}

	std::ofstream results;
	results.open(results_file);

	// Export transformation's parameters (i.e. n_min, beta_min/max, r)

	results << std::setprecision(std::numeric_limits<double>::max_digits10)
	           << "\nINPUTS:\n"
	           << "        lambda_min = "
	           << lambda_min
	           << "\n        lambda_max = "
	           << lambda_max;
	results << std::setprecision(std::numeric_limits<double>::max_digits10)
             << "\n\nOUTPUTS:\n"
             << "        number of nodes = "
             << n_min
             << "\n        beta_min = " 
             << beta_min
             << "\n        beta_max = "
             << beta_max
             << "\n        transformation order = "
             << r;

  // Export transformed G-L nodes

	std::ifstream nodes_file;
	nodes_file.open("../../data/TabulatedGlNodes.csv");

	std::string line_nodes, column_nodes;
	std::vector<float128> new_nodes, old_nodes;
	float128 transformed_node;
	int num_nodes;

	while(std::getline(nodes_file, line_nodes))
	{
		
		std::stringstream column_nodes_string(line_nodes);

		std::vector<std::string> row_nodes;

		while(std::getline(column_nodes_string, column_nodes, ','))
		{
			row_nodes.push_back(column_nodes);
		}

		num_nodes = stoi(row_nodes[0]);

		if(n_min == num_nodes)
		{

			results << "\n\nG-L NODES:\n"; 

			for(int k = num_nodes; k > num_nodes/2; k--)
			{

				old_nodes.push_back(static_cast<float128>(row_nodes[k]));
				transformed_node = boost::math::powm1(static_cast<float128>(row_nodes[k]), r) + 1;
				new_nodes.push_back(transformed_node);

        results << std::setprecision(std::numeric_limits<float128>::max_digits10)
                << "        "
                << -transformed_node
                << "\n";

			}

			for(int k = new_nodes.size() - 1; k >= 0; k--)
			{

				results << std::setprecision(std::numeric_limits<float128>::max_digits10)
                << "        "
                << new_nodes.at(k)
                << "\n";

			}

			break;

		}

	}

	nodes_file.close();

  // Export transformed G-L weights

	std::ifstream weights_file;
	weights_file.open("../../data/TabulatedGlWeights.csv");

	std::string line_weights, column_weights;
	std::vector<float128> new_weights;
	float128 transformed_weight;
	int num_weights;

	while(std::getline(weights_file, line_weights))
	{
		
		std::stringstream column_weights_string(line_weights);

		std::vector<std::string> row_weights;

		while(std::getline(column_weights_string, column_weights, ','))
		{
			row_weights.push_back(column_weights);
		}

		num_weights = stoi(row_weights[0].substr(1, row_weights[0].size() - 2));

		if(n_min == num_weights)
		{

			results << "\nG-L WEIGHTS:\n"; 

			for(int k = 1; k <= num_weights/2; k++)
			{

				transformed_weight = r*(boost::math::powm1(old_nodes.at(k-1), r-1) + 1)*static_cast<float128>(row_weights[k].substr(1, row_weights[k].size() - 2));
				new_weights.push_back(transformed_weight);

        results << std::setprecision(std::numeric_limits<float128>::max_digits10)
                << "        "
                << -transformed_weight
                << "\n";

			}

			for(int k = new_weights.size() - 1; k >= 0; k--)
			{

				results << std::setprecision(std::numeric_limits<float128>::max_digits10)
                << "        "
                << new_weights.at(k)
                << "\n";

			}

			break;

		}

	}

	weights_file.close();

	// Generate and export the plots if requested

	results.close();

	if(plot_flag.compare("y") == 0)
	{
		plots(n_min, beta_min, beta_max);
	}

}

void plots(const int& num_nodes, const double& min_exponent, const double& max_exponent)
{

	int print_startpoint = static_cast<int>(ceil(min_exponent)) -10;
	int print_endopoint = static_cast<int>(ceil(max_exponent)) + 10;
	double dp_epsilon = std::numeric_limits<double>::epsilon();

	std::ofstream epsilon;

	std::string mkdir = "mkdir -p ";
	std::string output_dir = "../../output/n=";
	std::string plot_subdir = mkdir + output_dir + std::to_string(num_nodes) + "/plot";

	system(plot_subdir.c_str());

	std::string epsilon_file = output_dir + std::to_string(num_nodes) + "/plot/Epsilon.csv";
	epsilon.open(epsilon_file.c_str());

	epsilon << std::setprecision(std::numeric_limits<float1k>::max_digits10) << print_startpoint << "," << dp_epsilon << "\n" << print_endopoint << "," << dp_epsilon;

	epsilon.close();

	std::ofstream f_exact, f_exact_env;

	std::string exact_file = output_dir  + std::to_string(num_nodes) + "/plot/ExactError.csv";
	f_exact.open(exact_file.c_str());
	std::string exact_env_file = output_dir  + std::to_string(num_nodes) + "/plot/ExactErrorEnv.csv";
	f_exact_env.open(exact_env_file.c_str());

	int n_samples = (max_exponent - min_exponent)*41;
	std::vector<double> lambda = linspace(static_cast<double>(print_startpoint), static_cast<double>(print_endopoint), n_samples);

	std::string envelope = "yes";
	std::string no_envelope = "no";

	std::cout << "\nComputing the error estimates with " << num_nodes << " quadrature nodes:	" << std::endl;

	for(int k = 0; k < n_samples - 2; k++)
	{
		
		f_exact << std::setprecision(std::numeric_limits<float1k>::max_digits10) << lambda[k+1] << "," << compute_estimate(lambda[k+1], num_nodes, no_envelope) << "\n";
		f_exact_env << std::setprecision(std::numeric_limits<float1k>::max_digits10) << lambda[k+1] << "," << compute_estimate(lambda[k+1], num_nodes, envelope) << "\n";

		print_progress(k, n_samples - 2);

	}

	std::string plot_dir = "bash ../../utilities/plot.sh ";
	std::string plot_input = plot_dir + std::to_string(num_nodes) + " " + std::to_string(print_startpoint) + " " + std::to_string(print_endopoint);

	system(plot_input.c_str());

	epsilon.close();
	f_exact.close();
	f_exact_env.close();

	std::cout << "\nComputation completed!" << std::endl;

}

void print_progress(const int& iter, const int& num_iter)
{

	int bar_length = 70;

	float ratio = static_cast<float>(iter)/static_cast<float>(num_iter);
	int progress = ceil(ratio*70);
	int progress_prctg = ceil(ratio*100.0);

	std::cout << "[";
	
	for(int k = 0; k <= bar_length; k++)
	{
		if(k <= progress)
		{
			std::cout << "#";
		}
		else
		{
			std::cout << ".";
		}
	}

	std::cout << "] " << std::setprecision(2) << progress_prctg << " %\r";
  std::cout.flush();

}