#include "gaussian_samples.h"
#include <math.h>
#include <chrono>
#include <gsl/gsl_integration.h>

// clear && rm outputs.txt && g++ -w curved_simplex.cpp -o curved_simplex -I../include -L. -lmtqr -lm -lquadmath -lgsl -lgslcblas && ./curved_simplex

template <typename precision>
float128 f_eval(int basis, precision x, precision y){
	// Evaluate the basis function in x and y
	float128 f_value;
	switch(basis){
		case 1:
		f_value = (-(pow(-1 + sqrt(1 - x),2)*(-1 - 4*(-1 + sqrt(2))*y + (-29 + 20*sqrt(2))*pow(y,2) + pow(x,2)*(-1 - 8*(-1 + sqrt(2))*y + 32*(-3 + 2*sqrt(2))*pow(y,2)) - 2*x*(-1 - 6*(-1 + sqrt(2))*y + 4*(-13 + 9*sqrt(2))*pow(y,2))))/(4 + 8*(-1 + sqrt(2))*x + 8*(-1 + sqrt(2))*y));
		break;

		case 2:
		f_value = ((-1 + sqrt(1 - x))*(2*(-1 + sqrt(1 - x)) + (5 - 6*sqrt(2) + 4*sqrt(2 - 2*x) - 2*sqrt(1 - x))*y + (-55 + 38*sqrt(2) - 36*sqrt(2 - 2*x) + 52*sqrt(1 - x))*pow(y,2) + 2*(-1 + 2*sqrt(1 - x))*pow(x,3)*(1 - sqrt(2) + 8*(-3 + 2*sqrt(2))*y) + x*(2 + 2*sqrt(2) - 4*sqrt(2 - 2*x) + (11 - 2*sqrt(2) + 28*sqrt(2 - 2*x) - 46*sqrt(1 - x))*y + 2*(89 - 61*sqrt(2) + 50*sqrt(2 - 2*x) - 74*sqrt(1 - x))*pow(y,2)) - 2*pow(x,2)*(-1 + 2*sqrt(2) - 4*sqrt(2 - 2*x) + 3*sqrt(1 - x) + 4*(8 - 5*sqrt(2) + 12*sqrt(2 - 2*x) - 18*sqrt(1 - x))*y + 8*(-3 + 2*sqrt(2))*(-3 + 2*sqrt(1 - x))*pow(y,2))))/(4.*pow(1 - x,1.5)*(1 + 2*(-1 + sqrt(2))*x + 2*(-1 + sqrt(2))*y));
		break;

		case 3:
		f_value = -((-1 + sqrt(1 - x))*(2 - 2*sqrt(1 - x) - 8*(-1 + sqrt(2))*(-1 + sqrt(1 - x))*y + (55 - 38*sqrt(2) + 36*sqrt(2 - 2*x) - 52*sqrt(1 - x))*pow(y,2) - 2*x*(2 - 2*sqrt(1 - x) - (-1 + sqrt(2))*(-11 + 10*sqrt(1 - x))*y + (89 - 61*sqrt(2) + 50*sqrt(2 - 2*x) - 74*sqrt(1 - x))*pow(y,2)) + 2*pow(x,2)*(1 - sqrt(1 - x) - (-1 + sqrt(2))*(-7 + 6*sqrt(1 - x))*y + 8*(-3 + 2*sqrt(2))*(-3 + 2*sqrt(1 - x))*pow(y,2))))/(4.*pow(1 - x,1.5)*(1 + 2*(-1 + sqrt(2))*x + 2*(-1 + sqrt(2))*y));
		break;

		case 4:
		f_value = (pow(2*(-1 + sqrt(2))*x*(-1 + x + y) + 2*(-1 + sqrt(1 - x))*(-1 + x)*(-1 + 2*(-1 + sqrt(2))*x - 2*(-1 + sqrt(2))*y),2) + pow((1 + 2*(-1 + sqrt(2))*x)*(-1 + x + y) + 2*(-1 + sqrt(1 - x))*(-1 + x)*(1 + 2*(-1 + sqrt(2))*x - 2*(-1 + sqrt(2))*y),2))/(4.*pow(1 - x,3)*(1 + 2*(-1 + sqrt(2))*x + 2*(-1 + sqrt(2))*y));
		break;

		case 5:
		f_value = -((-1 + 4*(-1 + sqrt(2))*(-1 + sqrt(1 - x))*(-1 + x) - 2*(-1 + sqrt(2))*x)*y*((1 + 2*(-1 + sqrt(2))*x)*(-1 + x + y) + 2*(-1 + sqrt(1 - x))*(-1 + x)*(1 + 2*(-1 + sqrt(2))*x - 2*(-1 + sqrt(2))*y)) + (2*(-1 + sqrt(2))*x*(-1 + x + y) + 2*(-1 + sqrt(1 - x))*(-1 + x)*(-1 + 2*(-1 + sqrt(2))*x - 2*(-1 + sqrt(2))*y))*(-2*(-1 + sqrt(2))*x*y + 2*(-1 + sqrt(1 - x))*(-1 + x)*(1 + 2*(-1 + sqrt(2))*y)))/(4.*pow(1 - x,3)*(1 + 2*(-1 + sqrt(2))*x + 2*(-1 + sqrt(2))*y));
		break;

		case 6:
		f_value = (pow(1 - 4*(-1 + sqrt(2))*(-1 + sqrt(1 - x))*(-1 + x) + 2*(-1 + sqrt(2))*x,2)*pow(y,2) + 4*pow((-1 + sqrt(2))*x*y - (-1 + sqrt(1 - x))*(-1 + x)*(1 + 2*(-1 + sqrt(2))*y),2))/(4.*pow(1 - x,3)*(1 + 2*(-1 + sqrt(2))*x + 2*(-1 + sqrt(2))*y));
		break;
	}

	// std::cout << std::setprecision(std::numeric_limits<float128>::digits10) << "f(x,y)=" << f_value << " @ x=" << x << ", y=" << y << "\n";

	return f_value;
}

class IntegrationWorkspace{
	// Wrapper for GSL's integration workspace
  	gsl_integration_workspace * wsp;

  	public:
  	IntegrationWorkspace(const size_t n=1000):
    	wsp(gsl_integration_workspace_alloc(n)) {}
  	~IntegrationWorkspace() { gsl_integration_workspace_free(wsp); }

  	void print_iter(){
  		std::cout << wsp -> size << "\n";
  	}

  	operator gsl_integration_workspace*() { return wsp; }
};

template <typename F>
class gsl_function_pp: public gsl_function {
	// Builds gsl_function from lambda
  	const F func;
  	static double invoke(double x, void *params) {
    return static_cast<gsl_function_pp*>(params)->func(x);
  	}
  	
  	public:
  	gsl_function_pp(const F& f) : func(f) {
    	function = &gsl_function_pp::invoke; //inherited from gsl_function
    	params   = this;                     //inherited from gsl_function
  	}
  	operator gsl_function*(){return this;}
};

template <typename F>
gsl_function_pp<F> make_gsl_function(const F& func){
	// Constructs integrand from template (helper function)
  	return gsl_function_pp<F>(func);
}

std::vector<double> integrate(){
	// Numerical cubature using QAGS
	double epsabs = 1e-13;
  	double epsrel = 1e-13;
  	size_t limit = 100;
  	double result, abserr, inner_result, inner_abserr;
  	std::vector<double> quadratures;

  	IntegrationWorkspace wsp1(limit);
  	IntegrationWorkspace wsp2(limit);

  	// BASIS 1
  	
	// Outer integral
	auto outer_1 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_1 = make_gsl_function([&](double y){
			return (-(pow(-1 + sqrt(1 - x),2)*(-1 - 4*(-1 + sqrt(2))*y + (-29 + 20*sqrt(2))*pow(y,2) + pow(x,2)*(-1 - 8*(-1 + sqrt(2))*y + 32*(-3 + 2*sqrt(2))*pow(y,2)) - 2*x*(-1 - 6*(-1 + sqrt(2))*y + 4*(-13 + 9*sqrt(2))*pow(y,2))))/(4 + 8*(-1 + sqrt(2))*x + 8*(-1 + sqrt(2))*y));}
		);
		gsl_integration_qags(inner_1, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_1, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Function basis 1:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quadratures.push_back(result);

	// BASIS 2
  	
	// Outer integral
	auto outer_2 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_2 = make_gsl_function([&](double y){
			return ((-1 + sqrt(1 - x))*(2*(-1 + sqrt(1 - x)) + (5 - 6*sqrt(2) + 4*sqrt(2 - 2*x) - 2*sqrt(1 - x))*y + (-55 + 38*sqrt(2) - 36*sqrt(2 - 2*x) + 52*sqrt(1 - x))*pow(y,2) + 2*(-1 + 2*sqrt(1 - x))*pow(x,3)*(1 - sqrt(2) + 8*(-3 + 2*sqrt(2))*y) + x*(2 + 2*sqrt(2) - 4*sqrt(2 - 2*x) + (11 - 2*sqrt(2) + 28*sqrt(2 - 2*x) - 46*sqrt(1 - x))*y + 2*(89 - 61*sqrt(2) + 50*sqrt(2 - 2*x) - 74*sqrt(1 - x))*pow(y,2)) - 2*pow(x,2)*(-1 + 2*sqrt(2) - 4*sqrt(2 - 2*x) + 3*sqrt(1 - x) + 4*(8 - 5*sqrt(2) + 12*sqrt(2 - 2*x) - 18*sqrt(1 - x))*y + 8*(-3 + 2*sqrt(2))*(-3 + 2*sqrt(1 - x))*pow(y,2))))/(4.*pow(1 - x,1.5)*(1 + 2*(-1 + sqrt(2))*x + 2*(-1 + sqrt(2))*y));}
  		);
		gsl_integration_qags(inner_2, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_2, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);
	
	std::cout << "Function basis 2:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";

	quadratures.push_back(result);

	// BASIS 3

	// Outer integral
	auto outer_3 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_3 = make_gsl_function( [&](double y){
			return -((-1 + sqrt(1 - x))*(2 - 2*sqrt(1 - x) - 8*(-1 + sqrt(2))*(-1 + sqrt(1 - x))*y + (55 - 38*sqrt(2) + 36*sqrt(2 - 2*x) - 52*sqrt(1 - x))*pow(y,2) - 2*x*(2 - 2*sqrt(1 - x) - (-1 + sqrt(2))*(-11 + 10*sqrt(1 - x))*y + (89 - 61*sqrt(2) + 50*sqrt(2 - 2*x) - 74*sqrt(1 - x))*pow(y,2)) + 2*pow(x,2)*(1 - sqrt(1 - x) - (-1 + sqrt(2))*(-7 + 6*sqrt(1 - x))*y + 8*(-3 + 2*sqrt(2))*(-3 + 2*sqrt(1 - x))*pow(y,2))))/(4.*pow(1 - x,1.5)*(1 + 2*(-1 + sqrt(2))*x + 2*(-1 + sqrt(2))*y));}
		);
		gsl_integration_qags(inner_3, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_3, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);
	
	std::cout << "Function basis 3:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";

	quadratures.push_back(result);

	// BASIS 4

	// Outer integral
	auto outer_4 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_4 = make_gsl_function( [&](double y){
			return (pow(2*(-1 + sqrt(2))*x*(-1 + x + y) + 2*(-1 + sqrt(1 - x))*(-1 + x)*(-1 + 2*(-1 + sqrt(2))*x - 2*(-1 + sqrt(2))*y),2) + pow((1 + 2*(-1 + sqrt(2))*x)*(-1 + x + y) + 2*(-1 + sqrt(1 - x))*(-1 + x)*(1 + 2*(-1 + sqrt(2))*x - 2*(-1 + sqrt(2))*y),2))/(4.*pow(1 - x,3)*(1 + 2*(-1 + sqrt(2))*x + 2*(-1 + sqrt(2))*y));}
		);
		gsl_integration_qags(inner_4, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_4, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);
	
	std::cout << "Function basis 4:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";

	quadratures.push_back(result);

	// BASIS 5

	// Outer integral
	auto outer_5 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_5 = make_gsl_function( [&](double y){
			return -((-1 + 4*(-1 + sqrt(2))*(-1 + sqrt(1 - x))*(-1 + x) - 2*(-1 + sqrt(2))*x)*y*((1 + 2*(-1 + sqrt(2))*x)*(-1 + x + y) + 2*(-1 + sqrt(1 - x))*(-1 + x)*(1 + 2*(-1 + sqrt(2))*x - 2*(-1 + sqrt(2))*y)) + (2*(-1 + sqrt(2))*x*(-1 + x + y) + 2*(-1 + sqrt(1 - x))*(-1 + x)*(-1 + 2*(-1 + sqrt(2))*x - 2*(-1 + sqrt(2))*y))*(-2*(-1 + sqrt(2))*x*y + 2*(-1 + sqrt(1 - x))*(-1 + x)*(1 + 2*(-1 + sqrt(2))*y)))/(4.*pow(1 - x,3)*(1 + 2*(-1 + sqrt(2))*x + 2*(-1 + sqrt(2))*y));}
		);
		gsl_integration_qags(inner_5, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_5, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);
	
	std::cout << "Function basis 5:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";

	quadratures.push_back(result);

	// BASIS 6

	// Outer integral
	auto outer_6 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_6 = make_gsl_function( [&](double y){
			return (pow(1 - 4*(-1 + sqrt(2))*(-1 + sqrt(1 - x))*(-1 + x) + 2*(-1 + sqrt(2))*x,2)*pow(y,2) + 4*pow((-1 + sqrt(2))*x*y - (-1 + sqrt(1 - x))*(-1 + x)*(1 + 2*(-1 + sqrt(2))*y),2))/(4.*pow(1 - x,3)*(1 + 2*(-1 + sqrt(2))*x + 2*(-1 + sqrt(2))*y));}
		);
		gsl_integration_qags(inner_6, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_6, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);
	
	std::cout << "Function basis 6:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";

	quadratures.push_back(result);

	// Return the results
  	return quadratures;
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
	int n_iter = 1;
	float128 one128 = 1.000000000000000000000000000000000;
	std::vector<std::vector<double>> outputs;
	std::vector<double> quadrature;

	// Outer quadrature samples (MTQR)
	auto start_mtqr = std::chrono::high_resolution_clock::now();
	for(int n=0; n<n_iter; n++){
		outputs = mtqr(0.0, 10.0);
	}
	auto stop_mtqr = std::chrono::high_resolution_clock::now();
	std::vector<double> outer_nodes = outputs[0];
	std::vector<double> outer_weights = outputs[1];

	// Inner quadrature samples (G-L)
	std::vector<std::vector<float128>> _outputs = remap(nodes[8], weights[8]);
	std::vector<float128> inner_nodes = _outputs[0];
	std::vector<float128> inner_weights = _outputs[1];

	// Exact integrals of the basis functions
	std::vector<double> integral = {0.0012241429923909649944, // basis function n°1
									0.0078405230975909855942,      // basis function n°2
									0.0084234049826403671523,      // basis function n°3
									0.17218405398906181066,        // basis function n°4
									0.10986262732342780534,        // basis function n°5
									0.13977718017066214501};       // basis function n°6

	// Numerically integrate the basis functions using GSL-QAGS
	auto start_qags = std::chrono::high_resolution_clock::now();
	for(int n=0; n<n_iter; n++){
		quadrature = integrate();
	}
	auto stop_qags = std::chrono::high_resolution_clock::now();

	// Print the time measurements
	auto t_mtqr = std::chrono::duration_cast<std::chrono::microseconds>(stop_mtqr - start_mtqr);
	auto t_qags = std::chrono::duration_cast<std::chrono::microseconds>(stop_qags - start_qags);
	std::cout << "Time of execution\n" 
			  << "  MTQR = " << t_mtqr.count()/n_iter << " 10e-6 seconds\n"
			  << "  QAGS = " << t_qags.count()/n_iter << " 10e-6 seconds\n";


	// Numerically integrate the basis functions using MTQR
	for(int b=0; b<integral.size(); b++){
		// Define the exact integral (using Mathematica's infinite symbolic precision)
		float128 exact = integral[b];
		
		float128 outer_quadrature = 0.0;
		// OUTER INTEGRAL (IN X BETWEEN 0 & 1)
		for(int n=0; n<outer_nodes.size(); n++){

			float128 inner_quadrature = 0.0;
			// INNER INTEGRAL (IN Y BETWEEN 0 & 1-X)
			for(int m=0; m<inner_nodes.size(); m++){
				float128 x = one128 - outer_nodes[n];
				if(outer_nodes[n] < 1e-32){
					x = one128 - 1e-32;
				}
				else{
					x = one128 - outer_nodes[n];
				}
				float128 y = static_cast<float128>(inner_nodes[m])*static_cast<float128>(outer_nodes[n]);
				inner_quadrature += inner_weights[m]*f_eval(b+1, x, y);
			}
			// Update the outer integral
			outer_quadrature += outer_weights[n]*outer_nodes[n]*inner_quadrature;
		}
		// Compare MTQR's and GSL-QAGS' results
		std::cout << std::setprecision(std::numeric_limits<double>::digits10)
											  << "Basis function " << b+1 << ": Mathematica's integral = " << exact
		                                      << "\n  MTQR vs Mathematica: " << abs(outer_quadrature-exact)/abs(exact)
											  << "\n  QAGS vs Mathematica: " << abs(quadrature[b]-exact)/abs(exact)
											  << "\n  MTQR vs QAGS: " << abs(outer_quadrature-quadrature[b])/abs(quadrature[b]) << "\n";
	}

	return 0;
}