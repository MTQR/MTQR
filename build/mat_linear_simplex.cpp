#include "gaussian_samples.h"
#include <math.h>
#include <chrono>
#include <gsl/gsl_integration.h>

// clear && rm outputs.txt && g++ -w mat_linear_simplex.cpp -o mat_linear_simplex -I../include -L. -lmtqr -lm -lquadmath -lgsl -lgslcblas && ./mat_linear_simplex

template <typename precision>
float128 f_eval(int row, int col, precision x, precision y){
	// Evaluate the basis function in x and y
	float128 f_value;
	switch(row){
		
		case 1:
		switch(col){
			
			case 1:
			f_value = pow(-1. + x,2) + pow(y,2);
			break;

			case 2:
			f_value = -x + pow(x,2) + (-1. + y)*y;
			break;

			case 3:
			f_value = (-1. + x)*x + pow(y,2);
			break;

			case 4:
			f_value = ((-1. + sqrt(1. - x))*(1. - 2*x + pow(x,2) + pow(y,2)))/2.;
			break;

			case 5:
			f_value = 1. - sqrt(1. - x) - x + y - y/sqrt(1. - x) - (y*(-1. + x + y))/(2.*pow(1. - x,1.5));
			break;

			case 6:
			f_value = 1. - sqrt(1. - x) - x - pow(y,2)/(2.*pow(1. - x,1.5));
			break;
		}
		break;

		case 2:
		switch(col){

			case 2:
			f_value = pow(x,2) + pow(-1. + y,2);
			break;

			case 3:
			f_value = pow(x,2) + (-1. + y)*y;
			break;

			case 4:
			f_value = ((-1. + sqrt(1. - x))*(-x + pow(x,2) + (-1. + y)*y))/2.;
			break;

			case 5:
			f_value = (-1. + 1./sqrt(1. - x))*x + (1. - y)*(-1. + 1./sqrt(1. - x) + (-1. + x + y)/(2.*pow(1. - x,1.5)));
			break;

			case 6:
			f_value = (-1. + 1./sqrt(1. - x))*x - ((-1. + y)*y)/(2.*pow(1. - x,1.5));
			break;
		}
		break;

		case 3:
		switch(col){

			case 3:
			f_value = pow(x,2) + pow(y,2);
			break;

			case 4:
			f_value = ((-1. + sqrt(1. - x))*(-x + pow(x,2) + pow(y,2)))/2.;
			break;

			case 5:
			f_value = (-1. + 1./sqrt(1. - x))*x + y - y/sqrt(1. - x) - (y*(-1. + x + y))/(2.*pow(1. - x,1.5));
			break;

			case 6:
			f_value = (-1. + 1./sqrt(1. - x))*x - pow(y,2)/(2.*pow(1. - x,1.5));
			break;
		}
		break;

		case 4:
		switch(col){

			case 4:
			f_value = (pow(-1. + sqrt(1. - x),2)*(1. - 2.*x + pow(x,2) + pow(y,2)))/4.;
			break;

			case 5:
			f_value = ((-1. + sqrt(1. - x))*(1. - sqrt(1. - x) - x + y - y/sqrt(1. - x) - (y*(-1. + x + y))/(2.*pow(1. - x,1.5))))/2.;
			break;

			case 6:
			f_value = ((-1. + sqrt(1. - x))*(-2.*(-1. + sqrt(1. - x) + x) - pow(y,2)/pow(1. - x,1.5)))/4.;
			break;
		}
		break;

		case 5:
		switch(col){

			case 5:
			f_value = pow(-1. + 1./sqrt(1. - x),2) + pow(-1. + 1./sqrt(1. - x) + (-1. + x + y)/(2.*pow(1. - x,1.5)),2);
			break;

			case 6:
			f_value = pow(-1. + 1./sqrt(1. - x),2) - (y*(1. - 2*sqrt(1. - x) + (-1. + 2*sqrt(1. - x))*x + y))/(4.*pow(-1. + x,3));
			break;
		}
		break;

		case 6:
		f_value = pow(-1. + 1./sqrt(1. - x),2) - pow(y,2)/(4.*pow(-1. + x,3));
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

std::vector<std::vector<double>> integrate(){
	// Numerical cubature using QAGS
	double epsabs = 1e-13;
  	double epsrel = 1e-13;
  	size_t limit = 100;
  	double result, abserr, inner_result, inner_abserr;
  	std::vector<std::vector<double>> quadratures;

  	IntegrationWorkspace wsp1(limit);
  	IntegrationWorkspace wsp2(limit);

  	// ROW N°1
  	std::vector<double> quad_1;

  	// ENTRY 1x1
  	
	// Outer integral
	auto outer_11 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_11 = make_gsl_function([&](double y){
			return pow(-1. + x,2) + pow(y,2);}
		);
		gsl_integration_qags(inner_11, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_11, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 1x1:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_1.push_back(result);

	// ENTRY 1x2
  	
	// Outer integral
	auto outer_12 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_12 = make_gsl_function([&](double y){
			return -x + pow(x,2) + (-1. + y)*y;}
		);
		gsl_integration_qags(inner_12, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_12, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 1x2:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_1.push_back(result);

	// ENTRY 1x3
  	
	// Outer integral
	auto outer_13 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_13 = make_gsl_function([&](double y){
			return (-1. + x)*x + pow(y,2);}
		);
		gsl_integration_qags(inner_13, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_13, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 1x3:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_1.push_back(result);

	// ENTRY 1x4
  	
	// Outer integral
	auto outer_14 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_14 = make_gsl_function([&](double y){
			return ((-1. + sqrt(1. - x))*(1. - 2*x + pow(x,2) + pow(y,2)))/2.;}
		);
		gsl_integration_qags(inner_14, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_14, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 1x4:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_1.push_back(result);

	// ENTRY 1x5
  	
	// Outer integral
	auto outer_15 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_15 = make_gsl_function([&](double y){
			return 1. - sqrt(1. - x) - x + y - y/sqrt(1. - x) - (y*(-1. + x + y))/(2.*pow(1. - x,1.5));}
		);
		gsl_integration_qags(inner_15, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_15, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 1x5:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_1.push_back(result);

	// ENTRY 1x6
  	
	// Outer integral
	auto outer_16 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_16 = make_gsl_function([&](double y){
			return 1. - sqrt(1. - x) - x - pow(y,2)/(2.*pow(1. - x,1.5));}
		);
		gsl_integration_qags(inner_16, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_16, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 1x6:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_1.push_back(result);
	quadratures.push_back(quad_1);

	// ROW N°2
	std::vector<double> quad_2;

	// ENTRY 2x2
  	
	// Outer integral
	auto outer_22 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_22 = make_gsl_function([&](double y){
			return pow(x,2) + pow(-1. + y,2);}
		);
		gsl_integration_qags(inner_22, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_22, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 2x2:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_2.push_back(result);

	// ENTRY 2x3
  	
	// Outer integral
	auto outer_23 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_23 = make_gsl_function([&](double y){
			return pow(x,2) + (-1. + y)*y;}
		);
		gsl_integration_qags(inner_23, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_23, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 2x3:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_2.push_back(result);

	// ENTRY 2x4
  	
	// Outer integral
	auto outer_24 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_24 = make_gsl_function([&](double y){
			return ((-1. + sqrt(1. - x))*(-x + pow(x,2) + (-1. + y)*y))/2.;}
		);
		gsl_integration_qags(inner_24, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_24, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 2x4:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_2.push_back(result);

	// ENTRY 2x5
  	
	// Outer integral
	auto outer_25 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_25 = make_gsl_function([&](double y){
			return (-1. + 1./sqrt(1. - x))*x + (1. - y)*(-1. + 1./sqrt(1. - x) + (-1. + x + y)/(2.*pow(1. - x,1.5)));}
		);
		gsl_integration_qags(inner_25, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_25, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 2x5:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_2.push_back(result);

	// ENTRY 2x6
  	
	// Outer integral
	auto outer_26 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_26 = make_gsl_function([&](double y){
			return (-1. + 1./sqrt(1. - x))*x - ((-1. + y)*y)/(2.*pow(1. - x,1.5));}
		);
		gsl_integration_qags(inner_26, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_26, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 2x6:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_2.push_back(result);
	quadratures.push_back(quad_2);

	// ROW N°3
	std::vector<double> quad_3;

	// ENTRY 3x3
  	
	// Outer integral
	auto outer_33 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_33 = make_gsl_function([&](double y){
			return pow(x,2) + pow(y,2);}
		);
		gsl_integration_qags(inner_33, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_33, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 3x3:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_3.push_back(result);

	// ENTRY 3x4
  	
	// Outer integral
	auto outer_34 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_34 = make_gsl_function([&](double y){
			return ((-1. + sqrt(1. - x))*(-x + pow(x,2) + pow(y,2)))/2.;}
		);
		gsl_integration_qags(inner_34, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_34, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 3x4:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_3.push_back(result);

	// ENTRY 3x5
  	
	// Outer integral
	auto outer_35 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_35 = make_gsl_function([&](double y){
			return (-1. + 1./sqrt(1. - x))*x + y - y/sqrt(1. - x) - (y*(-1. + x + y))/(2.*pow(1. - x,1.5));}
		);
		gsl_integration_qags(inner_35, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_35, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 3x5:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_3.push_back(result);

	// ENTRY 3x6
  	
	// Outer integral
	auto outer_36 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_36 = make_gsl_function([&](double y){
			return (-1. + 1./sqrt(1. - x))*x - pow(y,2)/(2.*pow(1. - x,1.5));}
		);
		gsl_integration_qags(inner_36, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_36, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 3x6:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_3.push_back(result);
	quadratures.push_back(quad_3);

	// ROW N° 4
	std::vector<double> quad_4;

	// ENTRY 4x4
  	
	// Outer integral
	auto outer_44 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_44 = make_gsl_function([&](double y){
			return (pow(-1. + sqrt(1. - x),2)*(1. - 2.*x + pow(x,2) + pow(y,2)))/4.;}
		);
		gsl_integration_qags(inner_44, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_44, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 4x4:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_4.push_back(result);

	// ENTRY 4x5
  	
	// Outer integral
	auto outer_45 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_45 = make_gsl_function([&](double y){
			return ((-1. + sqrt(1. - x))*(1. - sqrt(1. - x) - x + y - y/sqrt(1. - x) - (y*(-1. + x + y))/(2.*pow(1. - x,1.5))))/2.;}
		);
		gsl_integration_qags(inner_45, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_45, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 4x5:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_4.push_back(result);

	// ENTRY 4x6
  	
	// Outer integral
	auto outer_46 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_46 = make_gsl_function([&](double y){
			return ((-1. + sqrt(1. - x))*(-2.*(-1. + sqrt(1. - x) + x) - pow(y,2)/pow(1. - x,1.5)))/4.;}
		);
		gsl_integration_qags(inner_46, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_46, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 4x6:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_4.push_back(result);
	quadratures.push_back(quad_4);

	// ROW N°5
	std::vector<double> quad_5;

	// ENTRY 5x5
  	
	// Outer integral
	auto outer_55 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_55 = make_gsl_function([&](double y){
			return pow(-1. + 1./sqrt(1. - x),2) + pow(-1. + 1./sqrt(1. - x) + (-1. + x + y)/(2.*pow(1. - x,1.5)),2);}
		);
		gsl_integration_qags(inner_55, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_55, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 5x5:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_5.push_back(result);

	// ENTRY 5x6
  	
	// Outer integral
	auto outer_56 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_56 = make_gsl_function([&](double y){
			return pow(-1. + 1./sqrt(1. - x),2) - (y*(1. - 2*sqrt(1. - x) + (-1. + 2*sqrt(1. - x))*x + y))/(4.*pow(-1. + x,3));}
		);
		gsl_integration_qags(inner_56, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_56, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 5x6:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_5.push_back(result);
	quadratures.push_back(quad_5);

	// ROW N°6
	std::vector<double> quad_6;

	// ENTRY 6x6
  	
	// Outer integral
	auto outer_66 = make_gsl_function([&](double x){
		// Inner integral
		auto inner_66 = make_gsl_function([&](double y){
			return pow(-1. + 1./sqrt(1. - x),2) - pow(y,2)/(4.*pow(-1. + x,3));}
		);
		gsl_integration_qags(inner_66, 0, 1-x, epsabs, epsrel, limit, wsp1, &inner_result, &inner_abserr);
		return inner_result;
	});
	gsl_integration_qags(outer_66, 0, 1, epsabs, epsrel, limit, wsp2, &result, &abserr);

	std::cout << "Matrix entry 6x6:\n" << "  * inner iterations: ";
	wsp1.print_iter();
	std::cout << "  * outer iterations: ";
	wsp2.print_iter();
	std::cout << "\n";
	
	quad_6.push_back(result);
	quadratures.push_back(quad_6);

	// Return the compute quadratures
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
	std::vector<std::vector<double>> quadratures;

	// Outer quadrature samples (MTQR)
	auto start_mtqr = std::chrono::high_resolution_clock::now();
	for(int n=0; n<n_iter; n++){
		outputs = mtqr(0.0, 4.0);
	}
	auto stop_mtqr = std::chrono::high_resolution_clock::now();
	std::vector<double> outer_nodes = outputs[0];
	std::vector<double> outer_weights = outputs[1];

	// Inner quadrature samples (G-L)
	std::vector<std::vector<float128>> _outputs = remap(nodes[0], weights[0]);
	std::vector<float128> inner_nodes = _outputs[0];
	std::vector<float128> inner_weights = _outputs[1];

	// Exact integrals of the mass matrix
	std::vector<std::vector<float128>> integrals={// ROW N° 1
												 {0.33333333333333333333,   // matrix entry 1x1
												  -0.16666666666666666667,  // matrix entry 1x2
												  0.000000000000000000000,  // matrix entry 1x3
												  -0.018518518518518518519, // matrix entry 1x4
												  -0.066666666666666666667, // matrix entry 1x5
												  -0.13333333333333333333}, // matrix entry 1x6
												// ROW N° 2
												 {0.33333333333333333333,   // matrix entry 2x2
												  0.000000000000000000000,  // matrix entry 2x3
												  0.017195767195767195767,  // matrix entry 2x4
												  0.100000000000000000000,  // matrix entry 2x5
												  0.200000000000000000000}, // matrix entry 2x6
												// ROW N° 3
												 {0.16666666666666666667,   // matrix entry 3x3
												  0.0052910052910052910053, // matrix entry 3x4
												  0.100000000000000000000,  // matrix entry 3x5
												  0.033333333333333333333}, // matrix entry 3x6
												// ROW N° 4
												 {0.0018518518518518518519, // matrix entry 4x4
												  0.011507936507936507937,  // matrix entry 4x5
												  0.015079365079365079365}, // matrix entry 4x6
												// ROW N° 5
												 {0.250000000000000000000,  // matrix entry 5x5
												  0.20833333333333333333},  // matrix entry 5x6
												// ROW N° 6
												 {0.250000000000000000000}  // matrix entry 6x6
	};

	// Numerically compute the integrands using GSL-QAGS
	auto start_qags = std::chrono::high_resolution_clock::now();
	for(int n=0; n<n_iter; n++){
		quadratures = integrate();
	}
	auto stop_qags = std::chrono::high_resolution_clock::now();

	// Numerically compute the integrands using MTQR
	auto start_appl_mtqr = std::chrono::high_resolution_clock::now();
	for(int j=0; j<6; j++){
		// std::vector<float128> integral = integrals[j];
		int it=0;
		for(int k=j; k<6; k++){
			// Define the exact integral (using Mathematica's infinite symbolic precision)
			// float128 exact = integral[it];
			// Define the numerically compute integral using GSL's QAGS
			float128 quadrature = quadratures[j][it];
			// Initialise the numerically computed integral using MTQR
			float128 outer_quadrature = 0.0;
			// Define the samples of the outer quadrature
			if(j<3 && k<3){// Use G-L samples

				// OUTER INTEGRAL (IN X BETWEEN 0 & 1)
				for(int n=0; n<inner_nodes.size(); n++){
					float128 inner_quadrature = 0.0;
					// INNER INTEGRAL (IN Y BETWEEN 0 & 1-X)
					for(int m=0; m<inner_nodes.size(); m++){
						float128 x = one128 - inner_nodes[n];
						if(inner_nodes[n] < 1e-32){
							x = one128 - 1e-32;
						}
						else{
							x = one128 - inner_nodes[n];
						}
						float128 y = static_cast<float128>(inner_nodes[m])*static_cast<float128>(inner_nodes[n]);
						inner_quadrature += inner_weights[m]*f_eval(j+1, k+1, x, y);
					}
				// Update the outer integral
				outer_quadrature += inner_weights[n]*inner_nodes[n]*inner_quadrature;
				}
			}
			else{// Use MTQR samples
				
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
						inner_quadrature += inner_weights[m]*f_eval(j+1, k+1, x, y);
					}
				// Update the outer integral
				outer_quadrature += outer_weights[n]*outer_nodes[n]*inner_quadrature;
				}
			}
			// Compare MTQR's and GSL-QAGS' results
			/*
			if(abs(exact)>0){
				std::cout << "Matrix entry (j=" << j+1 << ",k=" << k+1 << ")\n  MTQR vs Mathematica: " << abs(outer_quadrature - exact)/abs(exact)
												  				   	   << "\n  QAGS vs Mathematica: " << abs(quadrature - exact)/abs(exact)
												  				   	   << "\n  MTQR vs QAGS: " << abs(outer_quadrature - quadrature)/abs(quadrature) << "\n";
			}
			else{
				std::cout << "Matrix entry (j=" << j+1 << ",k=" << k+1 << ")\n  MTQR vs Mathematica: " << abs(outer_quadrature)
												  				   	   << "\n  QAGS vs Mathematica: " << abs(quadrature)
												  				   	   << "\n  MTQR vs QAGS: " << abs(outer_quadrature - quadrature)/abs(quadrature) << "\n";
			}
			*/
			// Update counter
			it++;
		}
	}
	auto stop_appl_mtqr = std::chrono::high_resolution_clock::now();


	// Print the time measurements
	auto t_mtqr = std::chrono::duration_cast<std::chrono::microseconds>(stop_mtqr - start_mtqr);
	auto a_mtqr = std::chrono::duration_cast<std::chrono::microseconds>(stop_appl_mtqr - start_appl_mtqr);
	auto t_qags = std::chrono::duration_cast<std::chrono::microseconds>(stop_qags - start_qags);
	std::cout << "Time of execution\n" 
			  << "  MTQR (construction) = " << t_mtqr.count()/n_iter << " 10e-6 seconds\n"
			  << "  MTQR (application)  = " << a_mtqr.count()/21 << " 10e-6 seconds\n"
			  << "  QAGS = " << t_qags.count()/n_iter << " 10e-6 seconds\n";

	return 0;
}