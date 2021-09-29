//------------------------------------------------------------------------------------------------------------------
// File:      utilities/Plot.cpp
//
// Library:   QUASIMONT - QUAdrature of SIngular polynomials using MONomial Transformations:
//                        a C++ library for high precision quadrature of singular polynomials of non-integer degree
//
// Authors:   Davide Papapicco, Guido Lombardi, PhD
//
// Institute: Politecnico di Torino
//            C.so Duca degli Abruzzi, 24 - Torino (TO), Italia
//            Department of Electronics and Telecommunications (DET)
//            Electromagnetic modelling and applications Research Group
//------------------------------------------------------------------------------------------------------------------

#include "Quasimont.h"


//******************************************************************************************************************
//******************************************************************************************************************
//――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
//
//       FUNCTION: plot([n_min, {beta_min, beta_max}])
//                
//          INPUT: - [n_min, {beta_min, beta_max}] = output of function 'retrieveMonData'
//					
//		   OUTPUT: no outputs
//
//    DESCRIPTION: the exact asymptotic estimates are computed for each value of number of nodes. During
//				   the computation, or at the end of the monomial quadrature rule, the user has the choice of
//                 visualising the behaviour of the quadrature error estimate vs the exponents values by plotting
//                 the graph for the number of (quadrature) nodes computed by the library.
//
//     DISCLAIMER: due to its very nature this routine significantly slows down the performance of the library
//                 since it has to compute an exact estimate of the error for thousands of values of the exponent
//                 between -1 and lambda_max. For this reason a progress bar is displayed to help the user to have 
//                 a visual feedback on the status of computation (see lines 90~129).


void plot(const int& num_nodes, const double& min_exponent, const double& max_exponent)
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
	std::vector<float128> lambda = linspace(static_cast<float128>(print_startpoint), static_cast<float128>(print_endopoint), n_samples);

	std::string envelope = "yes";
	std::string no_envelope = "no";

	std::cout << "\nComputing the error estimates with " << num_nodes << " quadrature nodes:	" << std::endl;

	for(int k = 0; k < n_samples - 2; k++)
	{
		
		f_exact << std::setprecision(std::numeric_limits<float1k>::max_digits10) << lambda[k+1] << "," << computeEstimate(lambda[k+1], num_nodes, no_envelope) << "\n";
		f_exact_env << std::setprecision(std::numeric_limits<float1k>::max_digits10) << lambda[k+1] << "," << computeEstimate(lambda[k+1], num_nodes, envelope) << "\n";

		printProgressBar(k, n_samples - 2);

	}

	std::string plot_dir = "bash ../utilities/plot.sh ";
	std::string plot_input = plot_dir + std::to_string(num_nodes) + " " + std::to_string(print_startpoint) + " " + std::to_string(print_endopoint);

	system(plot_input.c_str());

	epsilon.close();
	f_exact.close();
	f_exact_env.close();

	std::cout << "\nComputation completed!" << std::endl;
}//――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――


//******************************************************************************************************************
//******************************************************************************************************************
//――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
//
//       FUNCTION: printProgressBar(iterator, number_of_iterations)
//                
//          INPUT: - iterator = value of the incremending integer index in the 'for' loop
//			       - number_of_iterations = end-value of the iterator in the 'for' loop
//					
//		   OUTPUT: no outputs
//
//    DESCRIPTION: in a traditional for(int k=0; k<num_k; k++) this function prints a progress bar on the terminal
//				   line to represent the status of completion of the loop.


void printProgressBar(const int& iter, const int& num_iter)
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
}//――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――