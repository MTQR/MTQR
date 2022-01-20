#include "Quasimont.h"

int main()
{
  // TEST ON STATIC LIBRARY THIRD-PARTY INTEGRATION
  std::vector<float128> coeff_sequence = {PI, 3};
  std::vector<float128> muntz_sequence = {-1/E, 1.0/2.0};
  quasimont(muntz_sequence, coeff_sequence);

  std::ifstream nodes_txt, weights_txt;
  std::vector<float128> nodes, weights;
  float128 loaded_value;
  // LOAD-IN COMPUTED NODES
  nodes_txt.open("output/Nodes.txt");
  while(nodes_txt >> loaded_value)
  {
    nodes.push_back(loaded_value);
  }
  nodes_txt.close();
  // LOAD-IN COMPUTED WEIGHTS
  weights_txt.open("output/Weights.txt");
  while(weights_txt >> loaded_value)
  {
    weights.push_back(loaded_value);
  }
  weights_txt.close();

  // COMPUTE QUADRATURE LOADED NODES AND WEIGHTS
  float128 In = 0;
  for(int k=0; k < muntz_sequence.size(); k++)
  {
    float128 In_mon = 0;
    for(int j=0; j<nodes.size(); j++)
    {
      In_mon += weights[j]*pow(static_cast<float128>(nodes[j]),muntz_sequence[k]);
    }
    In += coeff_sequence[k]*In_mon;
  }
  // PRINT COMPUTED QUADRATURE
  std::cout << "\n\n ** I_n(p(x)) = " << In << "  [with reloaded parameters in double precision] **" << std::endl;
  return 0;
}