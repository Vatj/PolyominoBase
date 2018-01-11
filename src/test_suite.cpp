#include "xorshift.hpp"
#include <chrono>
#include <iostream>


int main(int argc, char* argv[]) {
  
  uint32_t iterations = std::stoi(argv[1]);
  xorshift RNG_Xor;
  std::random_device rd;
  std::mt19937_64 RNG_MT (rd());
  std::minstd_rand nex(rd());
  std::ranlux48_base ran(rd());
  std::uniform_int_distribution<int> dis(0,199);
  std::uniform_real_distribution<double> dis2(0,1);
  auto begin = std::chrono::high_resolution_clock::now();
 
  for(uint32_t i = 0; i < iterations; ++i)
    {
      //RNG_Xor();
      dis(RNG_MT);
      dis2(RNG_MT);
    }
  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  std::cout <<"MT: "<< duration << "ns total, average : " << static_cast<double>(duration) / iterations << "ns." << std::endl;

  begin = std::chrono::high_resolution_clock::now();
 
  for(uint32_t i = 0; i < iterations; ++i)
    {
      dis(RNG_Xor);
      dis2(RNG_Xor);
    }
  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  std::cout <<"XORSHIFT: "<< duration << "ns total, average : " << static_cast<double>(duration) / iterations << "ns." << std::endl;

  begin = std::chrono::high_resolution_clock::now();
 
  for(uint32_t i = 0; i < iterations; ++i)
    {
      dis(nex);
      dis2(nex);
    }
  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  std::cout <<"Nex: "<< duration << "ns total, average : " << static_cast<double>(duration) / iterations << "ns." << std::endl;

  begin = std::chrono::high_resolution_clock::now();
 
  for(uint32_t i = 0; i < iterations; ++i)
    {
      dis(ran);
      dis2(ran);
    }
  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
  std::cout <<"Ran: "<< duration << "ns total, average : " << static_cast<double>(duration) / iterations << "ns." << std::endl;
}
