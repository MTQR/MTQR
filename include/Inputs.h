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

// lambda_max = (boost::math::constants::e<double>() + 0.25) = 2.968281828
// lambda_min = (-boost::math::constants::pi<double>()/4) = -0.785398163

std::cout << "\nInsert lambda_min: ";
std::cin >> lambda_min;

while(lambda_min <= -1)
{
    std::cout << "\n   ** ERROR ** lambda_min has to be (strictly) greater than -1 to be in a Muntz sequence\n\nInsert lambda_min: ";
    std::cin >> lambda_min;
}

std::cout << "Insert lambda_max: ";
std::cin >> lambda_max;

while(lambda_max <= lambda_min)
{
    std::cout << "\n   ** ERROR ** lambda_max has to be (strictly) greater than lambda_min\n\nInsert lambda_max: ";
    std::cin >> lambda_max;
}

std::string plot;

std::cout << "\nDo you want to plot the error estimates for your choice [y/n]? : ";
std::cin >> plot;