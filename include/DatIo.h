//------------------------------------------------------------------------------------------------------------------
// File:      include/DatIo.h
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

#ifndef DATIO_H
#define DATIO_H

// (SEE LINES 22~31 IN 'src/DatIo.cpp') Generate tabulated values for beta_min/beta_max for each even integer number nodes in [10,100]
void generateTabData();

// (SEE LINES 109~128 IN 'src/DatIo.cpp') Takes user-defined inputs from file
std::tuple<int, std::vector<float128>> getInputData(std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence);

// (SEE LINES 258~271 IN 'src/DatIo.cpp') Extract the values of beta_min and beta_max according to the computed minimum number of nodes
std::tuple<int, std::vector<double>> retrieveMonData(const int& comp_num_nodes);

// (SEE LINES 258~271 IN 'src/DatIo.cpp')
std::vector<double> collectData(const std::vector<double>& interval, const std::tuple<int, std::vector<float128>>& input_data, const double& transf_order);

// (SEE LINES 335~ IN 'src/DatIo.cpp') Degrade the precision of the new G-L nodes and weights to establish minimum data-type for double precision quadrature
void degradeData(std::tuple<double, std::vector<float50>, std::vector<float50>, std::vector<float50>, std::vector<float50>>& quad_params, std::vector<float128>& muntz_sequence, std::vector<float128>& coeff_sequence, const std::vector<double>& collected_data);

// (SEE LINES 335~349 IN 'src/DatIo.cpp') Computes and exports the resulting G-L weights and nodes aling with other ouputs
template<typename type>
void exportData(const std::tuple<double, std::vector<float50>, std::vector<float50>, std::vector<float50>, std::vector<float50>>& quad_params, const std::vector<type>& output_data, const std::vector<double>& collected_data);

#endif // DATIO_H