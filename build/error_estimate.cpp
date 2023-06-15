#include "mtqr.h"
#include <boost/multiprecision/mpfr.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/sin_pi.hpp>

// clear && rm -r error_estimate && g++ error_estimate.cpp -o error_estimate -I../include -L. -lmtqr -lm -lquadmath -lgsl -lgslcblas && ./error_estimate

typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<50>> float50;

using boost::math::policies::domain_error;
using boost::math::policies::errno_on_error;

typedef boost::math::policies::policy<
   boost::math::policies::domain_error<boost::math::policies::errno_on_error>,
   boost::math::policies::domain_error<boost::math::policies::errno_on_error>,
   boost::math::policies::domain_error<boost::math::policies::errno_on_error>,
   boost::math::policies::domain_error<boost::math::policies::errno_on_error>
> gamma_policy;

std::vector<float50> linspacedVector(const float50& start_type, const float50& end_type, const int& num_steps)
{
  std::vector<float50> linspaced;

  float50 start = static_cast<float50>(start_type);
  float50 end = static_cast<float50>(end_type);
  float50 steps = static_cast<float50>(num_steps);

  if (steps == 0)
  { 
    return linspaced;
  }

  if (steps == 1) 
  {
    linspaced.push_back(start);
    return linspaced;
  }

  float50 delta = (end - start) / (steps - 1);
  for(int k = 0; k < steps - 1; k++)
  {
    linspaced.push_back(start + delta*k);
  }
  linspaced.push_back(end);
  return linspaced;
}

template<typename T>
T aPrioriAsympEstimate(const T& lambda, const int& num_nodes, bool envelope)
{
  T common_factor = 2*(pow(static_cast<T>(2.0),-(static_cast<T>(1.0)+lambda)));
  T denominator = 2*num_nodes + lambda;

  T b1 = (boost::math::tgamma(static_cast<T>(2.0)*lambda, gamma_policy())*boost::math::tgamma(static_cast<T>(2.0*num_nodes) - lambda, gamma_policy()))/boost::math::tgamma(static_cast<T>(2.0)*lambda + static_cast<T>(2.0*num_nodes) - lambda, gamma_policy());
  T b2 = (boost::math::tgamma(static_cast<T>(2.0)*lambda, gamma_policy())*boost::math::tgamma(static_cast<T>(2.0) + static_cast<T>(2.0*num_nodes) - lambda, gamma_policy()))/boost::math::tgamma(static_cast<T>(2.0)*lambda + static_cast<T>(2.0) + static_cast<T>(2.0*num_nodes) - lambda, gamma_policy());

  T exact = common_factor*(pow(static_cast<T>(2.0), -lambda))*lambda*((b1/denominator) - b2/(2 + denominator));
  
  if(lambda + static_cast<T>(0.5) - 2*num_nodes >= 0 || envelope)
  {
    return fabs(exact*boost::math::sin_pi(lambda));
  }
  else
  {
    return fabs(exact);
  }
}
template float50 aPrioriAsympEstimate<float50>(const float50& input_lambda, const int& num_nodes, bool envelope);

template<typename T>
void plot(const int& num_nodes, const T& beta_min, const T& beta_max)
{
  std::cout << "\n\n *** Plotting the a-priori asymptotic error estimate with N=" << num_nodes << " ***" << std::endl;
  
  int print_startpoint = -1;
  int print_endopoint = static_cast<int>(ceil(beta_max)) + 10;

  std::ofstream epsilon;
  std::string epsilon_file = "./estimate/Epsilon.csv";
  epsilon.open(epsilon_file.c_str());
  epsilon << std::setprecision(std::numeric_limits<T>::max_digits10) << print_startpoint << "," << EPS << "\n" << print_endopoint << "," << EPS;
  epsilon.close();

  std::string csv = ".csv";
  std::string error_file = "./estimate/EstimateN=" + std::to_string(num_nodes) + csv;
  std::string env_error_file = "./estimate/EnvelopedEstimateN=" + std::to_string(num_nodes) + csv;
  std::ofstream error, env_error;
  error.open(error_file);
  env_error.open(env_error_file);

  int n_samples = (beta_max - beta_min)*100;
  std::vector<T> lambda = linspacedVector(static_cast<T>(print_startpoint), static_cast<T>(print_endopoint), n_samples);

  for(int k = 0; k < n_samples - 2; k++)
  {
    error << std::setprecision(std::numeric_limits<T>::max_digits10) << lambda[k+1] << "," << aPrioriAsympEstimate(lambda[k+1], num_nodes, true) << "\n";
    env_error << std::setprecision(std::numeric_limits<T>::max_digits10) << lambda[k+1] << "," << aPrioriAsympEstimate(lambda[k+1], num_nodes, false) << "\n";
    printProgressBar(k, n_samples - 2);
  }

  std::string plot_file = "bash ./plot.sh ";
  std::string cmd = plot_file + std::to_string(num_nodes) + " " + std::to_string(print_startpoint) + " " + std::to_string(print_endopoint);
  system(cmd.c_str());

  error.close();
  env_error.close();
}
template void plot<float50>(const int& num_nodes, const float50& beta_min, const float50& beta_max);

void printProgressBar(const int& iter, const int& num_iter)
{
  int bar_length = 70;

  float ratio = static_cast<float>(iter)/static_cast<float>(num_iter);
  int progress = ceil(ratio*70);
  int progress_prctg = ceil(ratio*100.0);

  std::cout << "[";
  
  for(int k = 0; k <= bar_length; k++)
  {
    if(k <= progress)
    {
      std::cout << "#";
    }
    else
    {
      std::cout << ".";
    }
  }

  std::cout << "] " << std::setprecision(2) << progress_prctg << " %\r";
  std::cout.flush();
}

int main(int argc, char** argv)
{
  int plot_nodes = 16;

  system("mkdir -p estimate");
  std::ofstream data("./estimate/TabBetas.csv");

  int num_samples = 15000, prev_j=0;
  std::vector<float50> lambda = linspacedVector(static_cast<float50>(-1.0), static_cast<float50>(1300.0), num_samples);

  float50 exact, tab_beta_min, tab_beta_max;
  int min_nodes = 5, max_nodes = 50;

  for(int n = min_nodes; n <= max_nodes; n++)
  {
    int nodes = 2*n;
    std::cout << "\nComputing for n = " << nodes << std::endl;

    for(int k = 0; k < num_samples - 2; k++)
    {
      std::cout << "    beta_min = " << lambda[k+1] << "\r";

      exact = aPrioriAsympEstimate(lambda[k+1], nodes, true);
      if(exact <= EPS)
      {
        std::cout << std::endl;
        tab_beta_min = lambda[k+1];

        if(prev_j < k+1)
        {
          prev_j = k+1;
        }
        for(int j = prev_j+1; j < num_samples - 2; j++)
        {
          std::cout << "    beta_max = " << lambda[j] << "\r";

          exact = aPrioriAsympEstimate(lambda[j], nodes, true);
          if(exact >= EPS)
          {
            tab_beta_max = lambda[j];
            prev_j = j;

            break;
          }
        }
        break;
      }
    }
    data << std::setprecision(std::numeric_limits<float50>::max_digits10) << nodes << "," << tab_beta_min << "," << tab_beta_max << "\n";

    if(nodes==plot_nodes)
    {
      plot(nodes, tab_beta_min, tab_beta_max);
      return 0;
    }
  }
  data.close();

  std::cout << "\n\n ** Data available in the 'estimate' subdirectory **\n ** Done! **\n" << std::endl;;
}