#ifndef UTILS_H
#define UTILS_H

// (SEE LINES 18~30 IN 'utilities/VecOps.cpp') Returns the float128 input vector in a type specified by the instantiation
template<typename type>
std::vector<type> castVector(const std::vector<float128>& input_vector, const type& type_infer);
// (SEE LINES 45~61 IN 'utilities/VecOps.cpp') Computes the inner product between two vectors (of the same type) avoiding numerical cancellation
template<typename type>
float128 doubleDotProduct(const std::vector<float128>& f_values, const std::vector<type>& weights);
// (SEE LINES 33~47 IN 'utilities/ErrTools.cpp') Generates n equispaced points between two input real numbers
template<typename type>
std::vector<type> linspacedVector(const type& start_type, const type& end_type, const int& num_steps);
// (SEE LINES 77~90 IN 'utilities/ErrTools.cpp')
template<typename type>
type aPrioriAsympEstimate(const type& input_lambda, const int& num_nodes);
// (SEE LINES 114~129 IN 'utilities/ErrTools.cpp')
template<typename type>
void plot(const int& num_nodes, const type& beta_min, const type& beta_max);
// (SEE LINES 167~180 IN 'utilities/ErrTools.cpp')
void printProgressBar(const int& iter, const int& num_iter);

#endif // UTILS_H