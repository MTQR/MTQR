//------------------------------------------------------------------------------------------------------------------
// File:      include/Quasimot.h
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

#ifndef QUASIMONT_H
#define QUASIMONT_H

#include <iostream>
#include <algorithm>
#include <iomanip>
#include <string>
#include <vector>
#include <tuple>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/powm1.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <gsl/gsl_poly.h>

typedef boost::multiprecision::float128 float128; // defines quadruple precision floating-point format datatype
typedef boost::multiprecision::mpf_float_50 float50; // defines 50 decimal digits floating-point format datatype
typedef boost::multiprecision::mpf_float_1000 float1k; // defines 1000 decimal digits floating-point format datatype

#include "Utils.h"  // includes header file for plotting and other utilities
#include "DatIo.h"  // includes header file for data I/O functions
#include "MonMap.h" // includes header file for functions computing the monomial map

void quasimont(std::vector<float128>& muntz_sequence, std::vector<float128>& poly_coeff, std::vector<double>& interval, std::string& plots);

#endif // QUASIMONT_H