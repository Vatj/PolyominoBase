#include <chrono>
#include <iostream>
#include <vector>
#include <random>
#include <fstream>

int main(int argc, char* argv[]) {
  
  std::random_device rd;
  std::mt19937 RNG_MT (rd());
  std::ofstream fout_size("Vectors", std::ios_base::out);
  for(int x=1;x<25;++x) {
    std::gamma_distribution<double> dis2(sqrt(1.*x),.5);
    for(int i=0;i<100000;++i)
      fout_size << dis2(RNG_MT)<<" ";
    std::cout<<dis2(RNG_MT)<<std::endl;
    fout_size << "\n";

  }

   
    
}
