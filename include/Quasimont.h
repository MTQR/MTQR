//---------------------------------------------------------------------------------------
// File:      include/Quasimont.h
//
// Library:   QUASIMONT-QUAdrature of SIngular polynomials using MONomial Transformations:
//                      a C++ library for high precision integration of singular 
//                      polynomials of non-integer degree
//
// Authors:   Guido Lombardi, Davide Papapicco
//
// Institute: Politecnico di Torino
//            C.so Duca degli Abruzzi, 24 - Torino (TO), Italia
//            Department of Electronics and Telecommunications (DET)
//            Electromagnetic modelling and applications Research Group
//---------------------------------------------------------------------------------------

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
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <gsl/gsl_poly.h>

namespace boomp = boost::multiprecision;

typedef boomp::float128 float128; // quadruple precision f.p. format
// typedef boomp::mpf_float_50 float50; // 50 decimal digits f.p. format
typedef boomp::number<boomp::mpfr_float_backend<50>> float50; // 50 decimal digits f.p. format

#include "Utils.h"  // includes header file for plotting and other utilities
#include "DatIo.h"  // includes header file for data I/O functions
#include "MonMap.h" // includes header file for functions computing the monomial map

#define EPS std::numeric_limits<double>::epsilon() // sets double machine-epsilon as treshold
#define PI boost::math::constants::pi<float50>() // defines pi with 50 decimal digits
#define E boost::math::constants::e<float50>() // defines e with 50 decimal digits

template<typename type>
void quasimont(std::vector<type>& muntz_sequence, std::vector<type>& coeff_sequence);

#endif // QUASIMONT_H