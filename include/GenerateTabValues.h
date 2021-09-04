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

std::cout << "\n\nComputation completed!\nTabulated data is now available in the 'NQL/monomial_transf/data' subdirectory.\n\nProgram terminated!\n" << std::endl;

return 0;