#ifndef _DATIO_H
#define _DATIO_H

// Generate tabulated values for beta_min/beta_max for each even integer number nodes in [10,100]
int generate_tab_values();

// Print the initial message (Result's taxonomy) and takes user-defined inputs from file
std::tuple<std::vector<double>, std::string> get_input_data(std::vector<double>& muntz_sequence);

// Extract the values of beta_min and beta_max according to the computed minimum number of nodes
std::tuple<int, std::vector<double>> extract_data(const int& comp_num_nodes);

// Computes and exports the resulting G-L weights and nodes aling with other ouputs
void export_results(const std::tuple<std::vector<double>, std::string>& user_data, const std::tuple<int, std::vector<double>>& data, const double& r);

// Plots the exact and enveloped error estimates given the number of quadrature nodes
void plots(const int& num_nodes, const double& min_exponent, const double& max_exponent);

// Prints a progress bar on-screen for the completion of the estimate computation
void print_progress(const int& iter, const int& num_iter);

#endif // _DATIO_H