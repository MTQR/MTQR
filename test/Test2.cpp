//------------------------------------------------------------------------------------------------------------------
// File:      test/test2.cpp
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
//          TITLE: Test 2
//
//     POLYNOMIAL: 
//
//     MOTIVATION:


int main(int argc, char** argv)
{
	std::vector<float128> coeff_sequence = {boost::math::constants::e<float128>(),
									  5.0,
									  -1.0,
									  1.0,
									  10.0};
	std::vector<float128> muntz_sequence = {boost::math::constants::e<float128>() + 0.25,
										  -boost::math::constants::pi<float128>()/4,
										  -0.5,
										  0,
										  2};
    std::vector<double> interval = {0, 1};
    std::string plots = "n";

	float128 I = quasimont(muntz_sequence, coeff_sequence, interval, plots);

	return 0;
}//――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――