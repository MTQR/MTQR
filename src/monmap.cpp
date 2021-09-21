#include "hiperquad.h"

template<typename type>
int compute_num_nodes(const type& lambda_max, const type& lambda_min)
{

	double dp_epsilon = std::numeric_limits<double>::epsilon();

	float1k l_max = static_cast<float1k>(lambda_max), l_min = static_cast<float1k>(lambda_min); 
	constexpr double c0 = -0.0040693, c1 = 0.00041296, d0 = 7.8147, d2 = 0.10123;

	float1k c2 = d2 + d2*l_min;
	float1k c3 = d0 + d0*l_min + l_min - l_max;
	float1k c4 = boost::math::powm1(1+l_max,3) + 1;

	double coeff7 = static_cast<double>(c1*(boost::math::powm1(c2,3)+1));
	double coeff6 = static_cast<double>(c0*(boost::math::powm1(c2,3)+1));
	double coeff5 = static_cast<double>(3*c1*(boost::math::powm1(c2,2) +1)*c3);
	double coeff4 = static_cast<double>(3*c0*(boost::math::powm1(c2,2) + 1)*c3);
	double coeff3 = static_cast<double>(3*c1*c2*(boost::math::powm1(c3,2)+1));
	double coeff2 = static_cast<double>(3*c0*c2*(boost::math::powm1(c3,2)+1));
	double coeff1 = static_cast<double>(c1*(boost::math::powm1(c3,3)+1));
	double coeff0 = static_cast<double>(c1*(boost::math::powm1(c3,3)+1) - c4);

  double coeff[8] = {coeff0, coeff1, coeff2, coeff3, coeff4, coeff5, coeff6, coeff7};
  double n[14];

  gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc	(8);
  gsl_poly_complex_solve(coeff, 8, w, n);
  gsl_poly_complex_workspace_free(w);

  int n_min, flag = 0;

  for (int k = 0; k < 7; k++)
  {
    if(fabs(n[2*k+1]) < dp_epsilon)
    {
    	if(n[2*k] > 0 && flag==0)
    	{
    		n_min = ceil(n[2*k]);
    		flag = 1;
    	}
    }
  }

	return n_min;

}
template int compute_num_nodes<double>(const double& lambda_max, const double& lambda_min); // Template mock instantiation for non-inline function

double compute_order(const std::vector<double>& lambdas, const std::vector<double>& betas)
{

	double lambda_min = lambdas[0], lambda_max = lambdas[1];
	double beta_min = betas[0], beta_max = betas[1];

	double r_max = (1 + beta_max)/(1 + lambda_max), r_min = (1 + beta_min)/(1 + lambda_min);
	double transf_order = (r_min + r_max)/2;
	
  std::cout << std::left 
            << std::setw(11) << ""
            << std::setw(11) << " r = " + std::to_string(transf_order) + ".\n\n";

  return transf_order;

}

template<typename type>
float1k compute_estimate(const type& input_lambda, const int& num_nodes, const std::string& envelope)
{

	float1k lambda = static_cast<float1k>(input_lambda);

	float1k common_factor = 2*(boost::math::powm1(2,-(1+lambda))+1);
	float1k denominator = 2*num_nodes + lambda;

  float1k b1 = (boost::math::tgamma(2*lambda)*boost::math::tgamma(2*num_nodes - lambda))/boost::math::tgamma(2*lambda + 2*num_nodes - lambda);
  float1k b2 = (boost::math::tgamma(2*lambda)*boost::math::tgamma(2 + 2*num_nodes - lambda))/boost::math::tgamma(2*lambda + 2 + 2*num_nodes - lambda);

  float1k exact = common_factor*(boost::math::powm1(2,-lambda)+1)*lambda*((b1/denominator) - b2/(2 + denominator));
  
  if(envelope.compare("yes") == 0)
	{
		if(lambda + 1/2 - 2*num_nodes >= 0)
	  {
	  	return fabs(exact*boost::math::sin_pi(lambda));
	  }
	  else
	  {
	  	return fabs(exact);
	  }
	}
	else
	{
		return fabs(exact*boost::math::sin_pi(lambda));
	}

}
template float1k compute_estimate<double>(const double& input_lambda, const int& num_nodes, const std::string& envelope); // Template mock instantiation for non-inline function

template<typename type>
std::vector<type> linspace(const type& start_type, const type& end_type, const int& num_steps)
{

  std::vector<type> linspaced;

  type start = static_cast<type>(start_type);
  type end = static_cast<type>(end_type);
  type steps = static_cast<type>(num_steps);

  if (steps == 0)
  { 
  	return linspaced;
  }

  if (steps == 1) 
  {
    linspaced.push_back(start);
    return linspaced;
  }

  type delta = (end - start) / (steps - 1);

  for(int k = 0; k < steps - 1; k++)
  {
    linspaced.push_back(start + delta*k);
  }

  linspaced.push_back(end);
  return linspaced;

}
template std::vector<double> linspace<double>(const double& start_type, const double& end_type, const int& num_steps); // Template mock instantiation for non-inline function