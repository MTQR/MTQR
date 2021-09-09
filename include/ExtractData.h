int n_min;
double beta_min, beta_max;

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

	if(n_min >= n_nodes)
	{

		beta_min = stof(row[1]);
		beta_max = stof(row[2]);

		std::cout << "\n\n   Results:\n";

		std::cout << std::left 
                  << std::setw(11) << ""
                  << std::setw(11) << "n_min = " + row[0] + ".\n";
        std::cout << std::left 
                  << std::setw(11) << ""
                  << std::setw(11)
                  << std::setprecision(std::numeric_limits<double>::max_digits10)
                  << "beta_min = " + row[1] + ".\n";
        std::cout << std::left 
                  << std::setw(11) << ""
                  << std::setprecision(std::numeric_limits<double>::max_digits10)
                  << std::setw(11) << "beta_max = " + row[2] + ".\n";

		break;

	}

}