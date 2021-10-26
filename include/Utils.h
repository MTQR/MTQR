//---------------------------------------------------------------------------------------
// File:      include/Utils.h
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

#ifndef UTILS_H
#define UTILS_H

// (SEE LINES 22~31 IN 'utilities/Utils.cpp') Returns the float50 input vector in a type specified by the instatiantiation
template<typename inputType, typename outputType>
void castVector(const std::vector<inputType>& input_vector, std::vector<outputType>& output_vector);

// (SEE LINES 22~31 IN 'utilities/Utils.cpp') Computes the inner product between two vectors (of the same type) avoiding numerical cancellation
template<typename type>
type orderedInnerProduct(const std::vector<type>& f_values, const std::vector<type>& weights);

// (SEE LINES 22~31 IN 'utilities/Utils.cpp') Generates n equispaced points between two input real numbers
template<typename type>
std::vector<type> linspace(const type& start_type, const type& end_type, const int& num_steps);

// Plots the exact and enveloped error estimates given the number of quadrature nodes
void plot(const int& num_nodes);//, const double& min_exponent, const double& max_exponent);

void printProgressBar(const int& iter, const int& num_iter);

#endif // UTILS_H