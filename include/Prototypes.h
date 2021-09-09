#ifndef PROTOTYPES_H
#define PROTOTYPES_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <math.h>
#include <boost/multiprecision/float128.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/powm1.hpp>
#include <boost/math/constants/constants.hpp>
#include <gsl/gsl_poly.h>

typedef boost::multiprecision::float128 float128; // defines quadruple precision floating-point format datatype
typedef boost::multiprecision::mpf_float_1000 float1k; // defines higher-than-quadruple precision (1000 decimal digits) floating-point format datatype

// Computes the number of minimum quadrature nodes by finding the real root of the 7-th degree polynomial equation in (62)
template<typename type>
int compute_n(const type& lambda_max, const type& lambda_min);

// Computes the exact estimate of the G-L quadrature error
template<typename type>
float1k compute_estimate(const type& input_lambda, const int& num_nodes, const std::string& envelope);

// Generates linearly spaced vector with n steps
template<typename type>
std::vector<type> linspace(const type& start_type, const type& end_type, const int& num_steps);

// Computes and exports the resulting G-L weights and nodes aling with other ouputs
void export_results(const double& lambda_min, const double& lambda_max, const int& n_min, const double& beta_min, const double& beta_max, const double& r, const std::string& plot_flag, const std::string& test_flag);

// Plots the exact and enveloped error estimates given the number of quadrature nodes
void plots(const int& num_nodes, const double& min_exponent, const double& max_exponent);

// Prints a progress bar on-screen for the completion of the estimate computation
void print_progress(const int& iter, const int& num_iter);

#endif