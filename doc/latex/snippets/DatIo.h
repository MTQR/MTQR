#ifndef DATIO_H
#define DATIO_H

// (SEE LINES 18~42 IN 'src/DatIo.cpp') Takes user-defined inputs from file
template<typename type>
std::tuple<int, std::vector<float128>> manageData(std::vector<type>& muntz_sequence, std::vector<type>& coeff_sequence);
// (SEE LINES 211~229 IN 'src/DatIo.cpp') Extract the values of beta_min and beta_max according to the computed minimum number of nodes
std::tuple<int, std::vector<float128>> streamMonMapData(const int& comp_num_nodes);
// (SEE LINES 274~294 IN 'src/DatIo.cpp') Degrade the precision of the new nodes and weights to establish minimum data-type for double precision quadrature
template<typename type>
void optimiseData(std::tuple<std::vector<float128>, std::vector<float128>, std::vector<float128>, std::vector<float128>>& quad_params, std::vector<type>& muntz_sequence, std::vector<type>& coeff_sequence);
// (SEE LINES 392~417 IN 'src/DatIo.cpp') Computes and exports the transformed weights and nodes along with other ouputs
template<typename type>
void exportNewData(const std::vector<type>& nodes, const std::vector<type>& weights, const std::vector<float128>& output_data);

#endif // DATIO_H