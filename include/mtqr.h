//---------------------------------------------------------------------------------------
// File:      include/mtqr.h
//
// Library:   MTQR - Monomial Transformation Quadrature Rule:
//                   a C++ library for high-precision integration of 
//                   generalised polynomials of non-integer degree
//
// Authors:   Guido Lombardi, Davide Papapicco
//
// Institute: Politecnico di Torino
//            C.so Duca degli Abruzzi, 24 - Torino (TO), Italia
//            Department of Electronics and Telecommunications (DET)
//            Electromagnetic modelling and applications Research Group
//---------------------------------------------------------------------------------------

#ifndef MTQR_H
#define MTQR_H

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
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <gsl/gsl_poly.h>

namespace boomp = boost::multiprecision;
typedef boomp::cpp_bin_float_quad float128; // quadruple precision f.p. format

#include "vector_operations.h"  // includes header file for vector data structures and operations
#include "data_management.h"  // includes header file for data handling and management 
#include "monomial_transformation.h" // includes header file for functions computing the monomial transformation

#define EPS std::numeric_limits<double>::epsilon() // sets double machine-epsilon as treshold
#define PI boost::math::constants::pi<float128>() // defines pi with 34 decimal digits
#define E boost::math::constants::e<float128>() // defines e with 34 decimal digits

// Loud mode
template<typename T>
void mtqr(std::vector<T>& muntz_sequence, std::vector<T>& coeff_sequence);
// Silent mode
std::vector<std::vector<double>> mtqr(double lambda_min, double lambda_max);

#endif // MTQR_H