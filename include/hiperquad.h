#ifndef _HIPERQUAD_H
#define _HIPERQUAD_H

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
typedef boost::multiprecision::mpf_float_1000 float1k; // defines higher-than-quadruple precision (1000 decimal digits) floating-point format datatype

#include "datio.h" // includes header file for data I/O functions
#include "monmap.h" // includes header file for computing the monomial map

void hiperquad(std::vector<double>& muntz_sequence);

#endif // _HIPERQUAD_H