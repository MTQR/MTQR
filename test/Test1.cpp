//------------------------------------------------------------------------------------------------------------------
// File:      test/test1.cpp
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
//          TITLE: Test 1
//
//     POLYNOMIAL: 
//
//     MOTIVATION:


int main(int argc, char** argv)
{
	std::vector<float128> coeff_sequence = {-1, 1};       // coeff = 1
	std::vector<float128> muntz_sequence = {-0.5, -boost::math::constants::pi<float128>()/4}; // lambda max -boost::math::constants::pi<float128>()/4 = -0.785398163397448309615660845819875699, n_min = 14
    std::vector<double> interval = {0, 1};
    std::string plots = "n";

	quasimont(muntz_sequence, coeff_sequence, interval, plots);

	return 0;
}//――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――