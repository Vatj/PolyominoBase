#include "GP_mapping.hpp"
#include <iostream>

/*External wrappers for python integration */
void LoadExistingTable(std::ifstream& fin,PhenotypeTable* pt_it);

extern "C"
{
  void GetPhenotypesIDs(const char* file_path_c,const char* file_name_c, uint8_t n_genes, uint8_t colours, bool file_of_genotypes);
  void PreProcessGenotypes(const char* file_path_c, uint8_t n_genes, uint8_t colours,bool file_of_genotypes);
}






