#ifndef VEC_OPS_H
#define VEC_OPS_H

// Returns the float128 input vector in a type specified by the instantiation
template<typename type>
std::vector<type> castVector(const std::vector<float128>& input_vector, const type& type_infer);

// Computes the inner product between two vectors (of the same type) avoiding numerical cancellation
template<typename type>
float128 doubleDotProduct(const std::vector<float128>& f_values, const std::vector<type>& weights);

#endif // VEC_OPS_H