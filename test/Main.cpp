//---------------------------------------------------------------------------------------
// File:      test/Main.cpp
//
// Library:   QUASIMONT-QUAdrature of SIngular polynomials using MONomial Transformations:
//                      a C++ library for high precision integration of singular 
//                      polynomials of non-integer degree
//
// Authors:   Guido Lombardi, PhD, Davide Papapicco
//
// Institute: Politecnico di Torino
//            C.so Duca degli Abruzzi, 24 - Torino (TO), Italia
//            Department of Electronics and Telecommunications (DET)
//            Electromagnetic modelling and applications Research Group
//---------------------------------------------------------------------------------------

#include "Quasimont.h"

int main(int argc, char** argv)
{

	std::vector<double> interval = {0, 1};

	std::vector<float128> coeff_sequence_bench_3 = {1,1};
	std::vector<float128> muntz_sequence_bench_3 = {17,35};

	quasimont(muntz_sequence_bench_3, coeff_sequence_bench_3, interval);

	std::vector<float128> coeff_sequence_bench_1 = {boost::math::constants::e<float128>(),
									  				5.0,
									  				-1.0,
									  				1.0,
									  				10.0};
	std::vector<float128> muntz_sequence_bench_1 = {boost::math::constants::e<float128>() + 0.25,
										  			-boost::math::constants::pi<float128>()/4,
										  			-0.5,
										  			0,
										  			2};

	quasimont(muntz_sequence_bench_1, coeff_sequence_bench_1, interval);

	std::vector<float128> coeff_sequence_bench_2 = {1};
	std::vector<float128> muntz_sequence_bench_2 = {-boost::math::constants::e<float128>()/3};

	quasimont(muntz_sequence_bench_2, coeff_sequence_bench_2, interval);

	return 0;
}