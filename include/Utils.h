//------------------------------------------------------------------------------------------------------------------
// File:      include/Utils.h
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

#ifndef UTILS_H
#define UTILS_H

//
template<typename type>
void castVector(const std::vector<float50>& input_vector, std::vector<type>& output_vector);

//
float1k orderedInnerProduct(std::vector<float1k>& input_vector);

// Generates linearly spaced vector with n steps
template<typename type>
std::vector<type> linspace(const type& start_type, const type& end_type, const int& num_steps);

// Plots the exact and enveloped error estimates given the number of quadrature nodes
void plot(const int& num_nodes, const double& min_exponent, const double& max_exponent);

// Prints a progress bar on-screen for the completion of the estimate computation
void printProgressBar(const int& iter, const int& num_iter);

#endif // UTILS_H