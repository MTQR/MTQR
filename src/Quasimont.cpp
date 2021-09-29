//------------------------------------------------------------------------------------------------------------------
// File:      src/Quasimot.cpp
// 
// Library:   QUASIMONT - QUAdrature of SIngular polynomials using MONomial Transformations:
//                        a C++ library for high precision integration of singular polynomials of non-integer degree
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
//       FUNCTION: post_map_integral = quasimont(muntz_sequence, coeff_sequence, I)
//                
//		    INPUT: - muntz_sequence = sequence of real exponents of the user-input polynomial
//			       - coeff_sequence = sequence of real coefficients of the user-input polynomial
//			       - I = [a,b] = interval of integration of the user-input polynomial (a < b)
//
//		   OUTPUT: - post_map_integral = output of function 'computeQuadGl'
//
//    DESCRIPTION: this routine takes the name of the library itself because it essentially instantiates and
//                 executes the ordered sequence of functions in it and it's interfaced direclty with the user-inputs.


void quasimont(std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence, std::vector<double>& interval, std::string& plots)
{
	// CHECK IF TABULATED VALUES EXIST AND COMPUTES THEM IF THEY DON'T
	std::ifstream datafile;
	datafile.open("../data/TabulatedErrorValues.csv");

	if(!datafile)
	{
		generateTabData();
	}

	// PRINT INITIAL MESSAGE AND SELECTS USER INPUTS
	auto input_data = getInputData(muntz_sequence, coeff_sequence);

	// EXTRACT N_MIN, BETA_MIN AND BETA_MAX
	auto monomial_data = retrieveMonData(std::get<0>(input_data)); // n_min = std::get<0>(input_data)

	// COMPUTE THE MONOMIAL TRANSFORMATION ORDER
	double transf_order = computeOrder(std::get<1>(input_data), std::get<1>(monomial_data)); // lambdas = std::get<1>(input_data), betas = std::get<1>(monomial_data)

	// COMPUTE AND EXPORT THE NEW G-L NODES & WEIGHTS
	auto quad_params = computeParams(transf_order, std::get<0>(monomial_data), interval); // n_min = std::get<0>(monomial_data)

	// COLLECT ALL THE OUTPUT DATA NECESSARY
	std::vector<double> collected_data = collectData(interval, input_data, transf_order);

	// CONVERTS AND EXPORTS NEW NODES AND WEIGHTS IN THE MINIMUM POSSIBLE FLOATING-POINT FORMAT
	degradeData(quad_params, muntz_sequence, coeff_sequence, collected_data);

	std::cout << "\n\nPROGRAM TERMINATED! Results are available in the 'output' subdirectory.\nPress any Key to exit ";
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	std::cin.get();
}//――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――


/*



*/