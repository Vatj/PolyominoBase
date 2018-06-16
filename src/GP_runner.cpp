#include "GP_runner.hpp"

#include <iterator>
#include <functional>
#include <set>

void LoadExistingTable(std::ifstream& fin,PhenotypeTable* pt_it) {
  std::string str;
  while (std::getline(fin, str)) {
    std::stringstream iss(str);
    int number;
    std::vector<uint8_t> phenotype_line;
    while (iss>>number)
      phenotype_line.push_back(static_cast<uint8_t>(number));
    Phenotype phen;
    phen.dx=phenotype_line[2];
    phen.dy=phenotype_line[3];
    phen.tiling=std::vector<uint8_t>(phenotype_line.begin()+4,phenotype_line.end());
    pt_it->known_phenotypes[phenotype_line[0]].emplace_back(phen);
  }
}

std::vector<Phenotype_ID> GetPhenotypeIDs(Genotype& genotype, uint8_t k_builds, PhenotypeTable* pt_it) {
  std::vector<Phenotype_ID> pIDs;
  Clean_Genome(genotype,0,false);
  const uint8_t n_genes=genotype.size()/4;
  std::map<uint8_t,uint8_t> dups=DuplicateGenes(genotype);
  std::map<uint8_t,Phenotype_ID> seed_m;
  for(uint8_t seed=0;seed<genotype.size()/4;++seed)
    seed_m[seed]=Stochastic::Analyse_Genotype_Outcome(genotype,k_builds,pt_it,seed);
  uint8_t index=0;
  Phenotype_ID phen_id;
  for(uint8_t gene=0;gene<n_genes; ++gene) {
    pIDs.emplace_back(dups.count(gene) ? seed_m[dups[gene]] : phen_id=seed_m[index++]);
  }
  std::sort(pIDs.begin(),pIDs.end());
  return pIDs;
}

void GetPhenotypesIDs(const char* file_path_c,const char* file_name_c, uint8_t n_genes, uint8_t colours, bool file_of_genotypes) {
  std::string file_path(file_path_c),file_name(file_name_c),str,details="_N"+std::to_string(n_genes)+"_C"+std::to_string(colours)+".txt";
  std::ifstream file_in(file_path+file_name);
  std::ofstream gfout(file_path+"Genotype_Codes"+details, std::ios_base::out);
  std::ofstream pfout(file_path+"Phenotype_Table"+details, std::ios_base::out);

  PhenotypeTable pt;
  uint8_t k_builds=10;
  Genotype genotype(n_genes*4);
  
  while (std::getline(file_in, str)) {
    if(file_of_genotypes) {
      std::istringstream is( str );
      genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );
    }
    else 
      index_to_genotype(std::stoull(str,nullptr),genotype,n_genes,colours);
    for(auto phen_id : GetPhenotypeIDs(genotype,k_builds,&pt)) 
      gfout<<+phen_id.first<<" "<<+phen_id.second<<" ";
    gfout<<"\n";
  }
  pt.PrintTable(pfout);
}


  
/****************/
/***** MAIN *****/
/****************/
int main(int argc, char* argv[]) {


  PhenotypeTable pt;
  Stochastic::STERIC_FORBIDDEN=true;
  for(auto x: Stochastic::AssemblePlasticGenotype({0,0,0,1, 0,0,1,2, 0,3,0,3, 0,4,0,0},&pt))
    std::cout<<+x.first<<" "<<x.second<<std::endl;
  Stochastic::STERIC_FORBIDDEN=false;
  for(auto x: Stochastic::AssemblePlasticGenotype({0,0,0,1, 0,0,1,2, 0,3,0,3, 0,4,0,0},&pt))
    std::cout<<+x.first<<" "<<x.second<<std::endl;
  return 0;
  
  if(argc>2) {
    std::cout<<argv[1][0]<<std::endl;

  }
  Genotype genotype{0,0,0,1,0,0,0,2};
  int x=0;
  for(Genotype neighbour : genotype_neighbourhood(genotype,2,std::stoi(argv[1]))) {
    ++x;
    for(auto q:neighbour)
      std::cout<<+q<<" ";
    std::cout<<std::endl;
  }
  std::cout<<x<<std::endl;
  return 0;

  
  Genotype geno,nullg;
  GenotypeGenerator ggenerator(std::stoi(argv[1]),std::stoi(argv[2]));
  ggenerator.init();
  while((geno=ggenerator())!=nullg) {
    for(auto x : geno)
      std::cout<<+x<<" ";
    std::cout<<std::endl;
  }
  
  

  return 0;


}












void PreProcessGenotypes(const char* file_path_c, uint8_t n_genes, uint8_t colours,bool file_of_genotypes) {
  std::string str,file_path(file_path_c),file_details="_N"+std::to_string(n_genes)+"_C"+std::to_string(colours);
  std::ifstream fin(file_path+"SampledGenotypes"+file_details+"_Iso.txt");
  std::ofstream fout(file_path+"SampledGenotypes"+file_details+"_Processed.txt");
  
  PhenotypeTable pt;
  uint8_t k_builds=5;
  std::map<Phenotype_ID,std::vector<Genotype> > phen_sets;
  Genotype genotype;
  
  while (std::getline(fin, str)) {
    if(file_of_genotypes) {
      std::istringstream is( str );
      genotype.assign( std::istream_iterator<int>( is ), std::istream_iterator<int>() );
    }
    const Genotype g_write(genotype);
    phen_sets[GetPhenotypeIDs(genotype, k_builds, &pt).front()].emplace_back(g_write);
    
  }

  for(auto kv : phen_sets) {
    for(Genotype& gen : kv.second) {
      for(uint8_t base : gen) {
	fout<<+base<<" ";
      }
      fout<<"\n";
    }
    fout<<"x\n";
  }
}
