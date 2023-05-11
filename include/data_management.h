//---------------------------------------------------------------------------------------
// File:      include/data_management.h
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

#ifndef DATA_MGT_H
#define DATA_MGT_H

// Takes user-defined inputs from file
template<typename type>
std::tuple<int, std::vector<float128>> manageData(std::vector<type>& muntz_sequence, std::vector<type>& coeff_sequence);

// Extract the values of beta_min and beta_max according to the computed minimum number of nodes
std::tuple<int, std::vector<float128>, int> streamMonMapData(const int& comp_num_nodes);

// Degrade the precision of the new nodes and weights to establish minimum data-type for double precision quadrature
template<typename type>
void optimiseData(std::tuple<std::vector<float128>, std::vector<float128>, std::vector<float128>, std::vector<float128>>& quad_params, std::vector<type>& muntz_sequence, std::vector<type>& coeff_sequence);

// Computes and exports the transformed weights and nodes along with other ouputs
template<typename type>
void exportNewData(const std::vector<type>& nodes, const std::vector<type>& weights, const std::vector<float128>& output_data);

#endif // DATA_MGT_H