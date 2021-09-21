#include "hiperquad.h"

void hiperquad(std::vector<double>& muntz_sequence)
{

	// CHECK IF TABULATED VALUES EXIST AND COMPUTES THEM IF THEY DON'T

	std::ifstream datafile;
	std::string tab_data = "../data/TabulatedErrorValues.csv";
	datafile.open(tab_data);

	if(!datafile)
	{
		int output = generate_tab_values();
	}

	// PRINT INITIAL MESSAGE AND SELECTS USER INPUTS

	auto user_data = get_input_data(muntz_sequence);
	std::vector<double> lambdas = std::get<0>(user_data); 
	double lambda_min = lambdas[0], lambda_max = lambdas[1];
	std::string plot = std::get<1>(user_data);

	// DERIVE MINIMUN NUMBER OF NODES

	int n_nodes = compute_num_nodes(lambda_max, lambda_min);

	// EXTRACT BETA_MIN AND BETA_MAX

	auto data = extract_data(n_nodes);
	int n_min = std::get<0>(data);
	std::vector<double> betas = std::get<1>(data);
	double beta_min = betas[0], beta_max = betas[1];

	// COMPUTE THE MONOMIAL TRANSFORMATION ORDER

	double transf_order = compute_order(lambdas, betas);

	// COMPUTE AND EXPORT THE NEW G-L NODES & WEIGHTS

	export_results(user_data, data, transf_order);

	std::cout << "\nPROGRAM TERMINATED!" << std::endl;

}