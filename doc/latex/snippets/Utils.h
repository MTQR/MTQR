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

// (SEE LINES 22~31 IN 'utilities/VecOps.cpp') Returns the float128 input vector in a type specified by the instatiantiation
template<typename type>
std::vector<type> castVector(const std::vector<float128>& input_vector, const type& type_infer);

// (SEE LINES 22~31 IN 'utilities/VecOps.cpp') Computes the inner product between two vectors (of the same type) avoiding numerical cancellation
template<typename type>
float128 doubleDotProduct(const std::vector<float128>& f_values, const std::vector<type>& weights);

// (SEE LINES 22~31 IN 'utilities/VecOps.cpp') Generates n equispaced points between two input real numbers
template<typename type>
std::vector<type> linspacedVector(const type& start_type, const type& end_type, const int& num_steps);

// (SEE LINES 22~41 IN 'utilities/ErrTools.cpp')
template<typename type>
type aPrioriAsympEstimate(const type& input_lambda, const int& num_nodes);

// (SEE LINES 22~41 IN 'utilities/ErrTools.cpp')
template<typename type>
void plot(const int& num_nodes, const type& beta_min, const type& beta_max);

// (SEE LINES 22~41 IN 'utilities/ErrTools.cpp')
void printProgressBar(const int& iter, const int& num_iter);

#endif // UTILS_H