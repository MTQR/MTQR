#include "hiperquad.h"

int main(int argc, char** argv)
{
	std::string pause;

	// TEST 1
	std::vector<double> muntz_sequence = {boost::math::constants::e<double>() + 0.25,
										  -boost::math::constants::pi<double>()/4,
										  -0.5,
										  0.0,
										  2.0};

	hiperquad(muntz_sequence);

	std::cout << "\nPress any key to continue"; 
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	std::cin.get();
	std::cout << "\n";
	

	/* TEST 2
	std::vector<double> muntz_sequence = {boost::math::constants::e<double>() + 0.25,
										  -boost::math::constants::pi<double>()/4,
										  -0.5,
										  0.0,
										  2.0};

	hiperquad(muntz_sequence);

	std::cout << "\nPress any key to continue"; 
	std::cin.get();

	// TEST 3
	std::vector<double> muntz_sequence = {boost::math::constants::e<double>() + 0.25,
										  -boost::math::constants::pi<double>()/4,
										  -0.5,
										  0.0,
										  2.0};

	hiperquad(muntz_sequence); */

	return 0;

}