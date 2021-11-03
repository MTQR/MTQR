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

// (SEE LINES 22~31 IN 'src/DatIo.cpp') Generates tabulated values for beta_min/beta_max for each even integer number nodes in [10,100]
void checkTabData();

// (SEE LINES 109~128 IN 'src/DatIo.cpp') Takes user-defined inputs from file
template<typename type>
std::tuple<int, std::vector<float50>> manageData(std::vector<type>& muntz_sequence, std::vector<type>& coeff_sequence);

// (SEE LINES 258~271 IN 'src/DatIo.cpp') Extract the values of beta_min and beta_max according to the computed minimum number of nodes
std::tuple<int, std::vector<float50>> streamMonMapData(const int& comp_num_nodes);

// (SEE LINES 335~ IN 'src/DatIo.cpp') Degrade the precision of the new G-L nodes and weights to establish minimum data-type for double precision quadrature
template<typename type>
void optimiseData(std::tuple<std::vector<float50>, std::vector<float50>, std::vector<float50>, std::vector<float50>>& quad_params, std::vector<type>& muntz_sequence, std::vector<type>& coeff_sequence);

// (SEE LINES 335~349 IN 'src/DatIo.cpp') Computes and exports the resulting G-L weights and nodes aling with other ouputs
template<typename type>
void exportNewData(const std::vector<type>& nodes, const std::vector<type>& weights, const std::vector<float50>& output_data);

#endif // DATIO_H