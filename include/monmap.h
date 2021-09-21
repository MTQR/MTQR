#ifndef _MONMAP_H
#define _MONMAP_H

// Computes the number of minimum quadrature nodes by finding the real root of the 7-th degree polynomial equation in (62)
template<typename type>
int compute_num_nodes(const type& lambda_max, const type& lambda_min);

// Computes the order (r) of the monomial map as a linear interpolation of r_min and r_max
double compute_order(const std::vector<double>& lambdas, const std::vector<double>& betas);

// Computes the exact estimate of the G-L quadrature error
template<typename type>
float1k compute_estimate(const type& input_lambda, const int& num_nodes, const std::string& envelope);

// Generates linearly spaced vector with n steps
template<typename type>
std::vector<type> linspace(const type& start_type, const type& end_type, const int& num_steps);

#endif // _MONMAP_H