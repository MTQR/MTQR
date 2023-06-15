#include "mtqr.h"

int main(){
  // Loud mode
  std::vector<float128> coeff_sequence = {3.0, PI};
  std::vector<float128> muntz_sequence = {1.0/E, -1.0/E};
  mtqr(muntz_sequence, coeff_sequence);

  std::ifstream nodes_txt, weights_txt;
  std::vector<float128> nodes, weights;
  float128 sample;
  
  nodes_txt.open("output/Nodes.txt");
  while(nodes_txt >> sample)
  {
    nodes.push_back(sample);
  }
  nodes_txt.close();
  
  weights_txt.open("output/Weights.txt");
  while(weights_txt >> sample)
  {
    weights.push_back(sample);
  }
  weights_txt.close();

  float128 In = 0.0;
  for(int k=0; k < muntz_sequence.size(); k++)
  {
    float128 In_mon = 0;
    for(int j=0; j<nodes.size(); j++)
    {
      In_mon += weights[j]*pow(static_cast<float128>(nodes[j]),muntz_sequence[k]);
    }
    In += coeff_sequence[k]*In_mon;
  }
  std::cout << "\n ** I_n(p(x)) = " << In << "  [reloaded samples in double precision]";
  
  // Silent mode
  double lambda_min = -1/double(E), lambda_max = 1/double(E);
  double coeff_min = double(PI), coeff_max = 3.0;
  std::vector<std::vector<double>> samples = mtqr(lambda_min, lambda_max);
  std::vector<double> dnodes = samples.at(0), dweights = samples.at(1);

  double dIn = 0.0;
  for(int j=0; j<dnodes.size(); j++){
  		dIn += dweights[j]*coeff_min*pow(dnodes[j],lambda_min)+coeff_max*pow(dnodes[j],lambda_max);
  	}
  std::cout << "\n ** I_n(p(x)) = " << In << "  [silent mode in double precision]\n\n";
  
  return 0;
}
