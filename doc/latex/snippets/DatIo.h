//---------------------------------------------------------------------------------------
// File:      include/DatIo.h
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

#ifndef DATIO_H
#define DATIO_H

// (SEE LINES 18~46 IN 'src/DatIo.cpp') Takes user-defined inputs from file
template<typename type>
std::tuple<int, std::vector<float128>> manageData(std::vector<type>& muntz_sequence, std::vector<type>& coeff_sequence);

// (SEE LINES 208~226 IN 'src/DatIo.cpp') Extract the values of beta_min and beta_max according to the computed minimum number of nodes
std::tuple<int, std::vector<float128>> streamMonMapData(const int& comp_num_nodes);

// (SEE LINES 268~292 IN 'src/DatIo.cpp') Degrade the precision of the new G-L nodes and weights to establish minimum data-type for double precision quadrature
template<typename type>
void optimiseData(std::tuple<std::vector<float128>, std::vector<float128>, std::vector<float128>, std::vector<float128>>& quad_params, std::vector<type>& muntz_sequence, std::vector<type>& coeff_sequence);

// (SEE LINES 372~399 IN 'src/DatIo.cpp') Computes and exports the resulting G-L weights and nodes aling with other ouputs
template<typename type>
void exportNewData(const std::vector<type>& nodes, const std::vector<type>& weights, const std::vector<float128>& output_data);

#endif // DATIO_H