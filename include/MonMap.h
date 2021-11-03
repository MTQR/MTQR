//---------------------------------------------------------------------------------------
// File:      include/MonMap.h
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

#ifndef MONMAP_H
#define MONMAP_H

// (SEE LINES 22~41 IN 'src/MonMap.cpp') Computes the exact estimate of the G-L quadrature error
float50 computeErrorEstimate(const float50& input_lambda, const int& num_nodes);

// (SEE LINES 79~93 IN 'src/MonMap.cpp') Computes the optimal lambda_max when the input polynomial is a monomial and a maximum number of nodes is required
float50 computeLambdaMax(float50& lambda_min, int num_nodes);

// (SEE LINES 114~129 IN 'src/MonMap.cpp') Computes the number of minimum quadrature nodes by finding the real root of the 7-th degree polynomial equation in (62)
int computeNumNodes(const float50& lambda_min, const float50& lambda_max);

// (SEE LINES 182~194 IN 'src/MonMap.cpp') Computes the order (r) of the monomial map as a linear interpolation of r_min and r_max
double computeMapOrder(const std::vector<float50>& lambdas, const std::vector<float50>& betas);

// (SEE LINES 218~234 IN 'src/MonMap.cpp') Computes the new nodes and weights of the G-L formula
std::tuple<std::vector<float50>, std::vector<float50>, std::vector<float50>, std::vector<float50>> computeParamsGl(const double& r, const int& n_min);

// (SEE LINES 349~365 IN 'src/MonMap.cpp') Computes the numerical integral with G-L quadrature formula
template<typename type1, typename type2>
float50 computeQuadGl(const std::vector<type1>& nodes, const std::vector<type1>& weights, std::vector<type2>& muntz_sequence, std::vector<type2>& coeff_sequence);

// (SEE LINES 450~462 IN 'src/MonMap.cpp') Computes the a-posteriori relative error of the G-L quadrature using the new nodes and weights
template<typename type>
float50 computeExactError(const float50& In, std::vector<type>& muntz_sequence, std::vector<type>& coeff_sequence);

#endif // MONMAP_H