#include "xorshift.hpp"
#include <chrono>
#include <iostream>


int main(int argc, char* argv[]) {
  
  uint32_t iterations = std::stoi(argv[1]);
  xorshift RNG_Xor;
  std::mt19937 RNG_MT (182751825);
  auto begin = std::chrono::high_resolution_clock::now();
 
  for(uint32_t i = 0; i < iterations; ++i)
    {
      //RNG_Xor();
      RNG_MT();
    }
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  std::cout << duration << "ns total, average : " << static_cast<double>(duration) / iterations << "ns." << std::endl;
 
}
