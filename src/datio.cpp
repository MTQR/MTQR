#include "hiperquad.h"

int generate_tab_values()
{

	std::ofstream data("../data/TabulatedErrorValues.csv");

	int num_samples = 15000, prev_j;
	double start = -1, end = 1300;
	std::vector<double> lambda = linspace(start, end, num_samples);

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

			exact = compute_estimate(lambda[k+1], nodes, envelope);

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

					exact = compute_estimate(lambda[j], nodes, envelope);

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

	return 0;

}

std::tuple<std::vector<double>, std::string> get_input_data(std::vector<double>& muntz_sequence)
{
	
	int w = 60;

	std::string dashed_line = "-";

	for(int k = 1; k < w; k++)
	{
	    dashed_line += "-";
	}

	std::cout << "//" + dashed_line + "\\\\" << std::endl;

	std::cout << "|" << std::right
	          << std::setw(w+4) << "|\n";
	std::cout << std::left 
	          << std::setw(ceil(w-21)/2) << "|"
	          << std::setw(ceil(w-21)/2) << "MONOMIAL TRANSORMATION RULE";
	std::cout << std::right
	          << std::setw(ceil(w-21)/2) << "|\n";

	std::cout << "|" << std::right
	          << std::setw(w+4) << "|\n";
	std::cout << "|" << std::right
	          << std::setw(w+4) << "|\n";

	std::cout << "|  Results taxonomy:";
	std::cout << std::right
	          << std::setw(w-15) << "|\n";
	std::cout << "|" << std::right
	          << std::setw(w+4) << "|\n";

	std::cout << std::left 
	          << std::setw(11) << "|"
	          << std::setw(11) << "n_min.............Minimum number of nodes.";
	std::cout << std::right
	          << std::setw(w-48) << "|\n";
	std::cout << std::left 
	          << std::setw(11) << "|"
	          << std::setw(11) << "beta_min..........Lower bound for the exponents.";
	std::cout << std::right
	          << std::setw(w-54) << "|\n";
	std::cout << std::left 
	          << std::setw(11) << "|"
	          << std::setw(11) << "beta_max..........Upper bound for the exponents.";
	std::cout << std::right
	          << std::setw(w-54) << "|\n";
	std::cout << std::left 
	          << std::setw(11) << "|"
	          << std::setw(11) << "r.................Transformation order.";
	std::cout << std::right
	          << std::setw(w-45) << "|\n";

	std::cout << "|" << std::right
	          << std::setw(w+4) << "|\n";
	std::cout << "|" << std::right
	          << std::setw(w+4) << "|\n";

	std::cout << "\\\\" + dashed_line + "//" << std::endl;

	double lambda_max, lambda_min;

	sort(muntz_sequence.begin(), muntz_sequence.end());

	if(muntz_sequence[0] <= -1.0)
	{
		std::cout << "\n   ** ERROR ** lambda_min has to be (strictly) greater than -1 to be in a Muntz sequence\n";
	}
	else
	{
		lambda_min = muntz_sequence[0]; 
	}

	lambda_max = muntz_sequence.back();

	std::cout << "\n ** Accepted sequence of Muntz exponents ** \n";
	std::cout << "    {" << muntz_sequence[0];
	for(int k=1; k < muntz_sequence.size(); k++)
	{
		std::cout << ", " << muntz_sequence[k];
	}
	std::cout << "}\n";

	std::vector<double> lambdas = {lambda_min, lambda_max};

	std::string plot;

	std::cout << "\n ** Do you want to plot the error estimates for your choice [y/n]? : ";
	std::cin >> plot;

	std::tuple<std::vector<double>, std::string> output = std::make_tuple(lambdas, plot);

	return output;

}

std::tuple<int, std::vector<double>> extract_data(const int& comp_num_nodes)
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

			std::cout << "\n\n    Results:\n";

			std::cout << std::left 
	                  << std::setw(11) << ""
	                  << std::setw(11) << " n_min = " + row[0] + ".\n";
	        std::cout << std::left 
	                  << std::setw(11) << ""
	                  << std::setw(11)
	                  << std::setprecision(std::numeric_limits<double>::max_digits10)
	                  << " beta_min = " + row[1] + ".\n";
	        std::cout << std::left 
	                  << std::setw(11) << ""
	                  << std::setprecision(std::numeric_limits<double>::max_digits10)
	                  << std::setw(11) << " beta_max = " + row[2] + ".\n";

			break;

		}

	}

	std::vector<double> betas = {beta_min, beta_max};

	return std::make_tuple(n_min, betas);

}

void export_results(const std::tuple<std::vector<double>, std::string>& user_data, const std::tuple<int, std::vector<double>>& data, const double& r)
{

	std::vector<double> lambdas = std::get<0>(user_data), betas = std::get<1>(data);
	double lambda_min = lambdas[0],lambda_max = lambdas[1], beta_min = betas[0], beta_max = betas[1];
	int n_min = std::get<0>(data);
	std::string plot = std::get<1>(user_data);

	std::string mkdir = "mkdir -p ";
	std::string output_dir = "output";
	std::string results_dir = mkdir + output_dir;
	system(results_dir.c_str());

	std::string results_file = output_dir + "/Results.txt";
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
	nodes_file.open("../data/TabulatedGlNodes.csv");

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
	weights_file.open("../data/TabulatedGlWeights.csv");

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

	if(plot.compare("y") == 0)
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

	std::string epsilon_file = "./output/Epsilon.csv";
	epsilon.open(epsilon_file.c_str());

	epsilon << std::setprecision(std::numeric_limits<float1k>::max_digits10) << print_startpoint << "," << dp_epsilon << "\n" << print_endopoint << "," << dp_epsilon;
	epsilon.close();

	std::ofstream f_exact, f_exact_env;

	f_exact.open("./output/ExactError.csv");
	f_exact_env.open("./output/ExactErrorEnv.csv");

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

	std::string plot_dir = "bash ../utilities/plot.sh ";
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