//---------------------------------------------------------------------------------------
// File:      include/vector_operations.h
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

#ifndef VEC_OPS_H
#define VEC_OPS_H

// Returns the float128 input vector in a type specified by the instantiation
template<typename type>
std::vector<type> castVector(const std::vector<float128>& input_vector, const type& type_infer);

// Computes the inner product between two vectors (of the same type) avoiding numerical cancellation
template<typename type>
float128 doubleDotProduct(const std::vector<float128>& f_values, const std::vector<type>& weights);

// Generates n equispaced points between two input real numbers
template<typename type>
std::vector<type> linspacedVector(const type& start_type, const type& end_type, const int& num_steps);

//
template<typename type>
type aPrioriAsympEstimate(const type& input_lambda, const int& num_nodes);

//
template<typename type>
void plot(const int& num_nodes, const type& beta_min, const type& beta_max);

//
void printProgressBar(const int& iter, const int& num_iter);

#endif // VEC_OPS_H