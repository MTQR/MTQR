#include "mtqr.h"
#include "curved_simplex.h"

// clear && g++ -w inverse_jacobian.cpp -o inverse_jacobian -I../include -L. -lmtqr -lm -lquadmath -lgsl -lgslcblas && ./inverse_jacobian

float128 f_eval(float128 x, float128 y){
	float128 f_value = 1/(1 + 2*(-1 + sqrt(2))*x + 2*(-1 + sqrt(2))*y);
	return f_value;
}

std::vector<std::vector<float128>> remap(std::vector<float128> nodes, std::vector<float128> weights){
	// Maps G-L nodes and weights from [-1,1] to [0,1]
	float128 a = 0.000000000000000000000000000000000;
	float128 b = 1.000000000000000000000000000000000;
	float128 jacobian = (b-a)/2;

	for(int j=0; j<nodes.size(); j++){
		nodes[j] = jacobian*nodes[j] + (a+b)/static_cast<float128>(2.000000000000000000000000000000000);
		weights[j] = jacobian*weights[j];
	}

	std::vector<std::vector<float128>> outputs = {nodes, weights};
	return outputs;
}

int main(int argc, char* argv[]){
	// Define the exact value of the integral
	float128 exact = {0.327806801937078208985509075433014};
	
	// Select different sets of G-L samples (from 2 to 18)
	for(int n=0; n<=16; n++){
		std::vector<std::vector<float128>> samples = remap(nodes[n], weights[n]);
		std::vector<float128> _nodes = samples[0];
		std::vector<float128> _weights = samples[1];

		// Outer quadrature in x between 0 and 1
		float128 outer = 0.000000000000000000000000000000000;
		for(int j=0; j<_nodes.size(); j++){
			// std::cout << "j=" << j << ": " << _nodes[j] << ", " << _weights[j] << "\n";
			float128 x = 1.000000000000000000000000000000000 - _nodes[j];
			// Inner quadrature in y between 0 and 1-x
			float128 inner = 0.000000000000000000000000000000000;
			for(int k=0; k<_nodes.size(); k++){
				// std::cout << "   k=" << k << ": " << _nodes[k] << ", " << _weights[k] << "\n";
				float128 y = _nodes[j]*_nodes[k];
				inner += _weights[k]*f_eval(x, y); 
			}
			outer += _weights[j]*_nodes[j]*inner;
		}

		// Print the results
		std::cout << std::setprecision(std::numeric_limits<float128>::digits10) << "  With " << n+2 << " samples:  " << abs(outer - exact)/abs(exact) << "\n";
	}
	return 0;
}