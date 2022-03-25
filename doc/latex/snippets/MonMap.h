#ifndef MONMAP_H
#define MONMAP_H

// (SEE LINES 18~37 IN 'src/MonMap.cpp') Computes the optimal lambda_max when the input polynomial is a monomial and a maximum number of nodes is required
float128 computeLambdaMax(float128& lambda_min, int num_nodes);
// (SEE LINES 53~70 IN 'src/MonMap.cpp') Computes the number of minimum quadrature nodes by finding the real root of the 7-th degree polynomial equation in (62)
int computeNumNodes(const float128& lambda_min, const float128& lambda_max);
// (SEE LINES 112~125 IN 'src/MonMap.cpp') Computes the order (r) of the monomial map as a linear interpolation of r_min and r_max
double computeMapOrder(const std::vector<float128>& lambdas, const std::vector<float128>& betas);
// (SEE LINES 184~202 IN 'src/MonMap.cpp') Computes the new nodes and weights of the monomial transformation quadrature rule
std::tuple<std::vector<float128>, std::vector<float128>, std::vector<float128>, std::vector<float128>> computeQuadParams(const double& r, const int& n_min);
// (SEE LINES 309~328 IN 'src/MonMap.cpp') Computes the numerical integral for a given quadrature rule
template<typename type1, typename type2>
float128 computeQuadRule(const std::vector<type1>& nodes, const std::vector<type1>& weights, std::vector<type2>& muntz_sequence, std::vector<type2>& coeff_sequence);
// (SEE LINES 356~378 IN 'src/MonMap.cpp') Computes the a-posteriori relative error of the quadrature rule
template<typename type>
float128 computeExactError(const float128& In, std::vector<type>& muntz_sequence, std::vector<type>& coeff_sequence, bool& print_primitive);

#endif // MONMAP_H