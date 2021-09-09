double lambda_max = boost::math::constants::pi<double>()*boost::math::constants::e<double>(); // 8.539734222673566
double lambda_min = -1.0/2.0;                                                                 // -0.5

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

w = 37;

std::string driver_line = "#";

for(int k = 1; k < w; k++)
{
    driver_line += "#";
}

std::cout << "\n";
std::cout << std::left
          << std::setw(14) << ""
          << std::setw(14) << driver_line
          << "\n";
std::cout << std::left
          << std::setw(14) << ""
          << std::setw(14) << "## Executing example driver: test3 ##"
          << "\n";
std::cout << std::left
          << std::setw(14) << ""
          << std::setw(14) << driver_line
          << "\n";

std::cout << "\nlambda_min = " << lambda_min;
std::cout << "\nlambda_max = " << lambda_max;

std::string plot = "n";
